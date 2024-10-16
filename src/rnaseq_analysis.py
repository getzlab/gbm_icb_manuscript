#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 09:46:38 2020
Modified in October 2021 by Conor Messer

@author: kschluet
"""

import sys
sys.path.insert(1, './src/')

import pickle
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import dalmatian
import re
import os
import math
import scipy.stats as ss
import itertools
from sklearn.decomposition import PCA
from IPython.display import Image
import matplotlib.pyplot as plt
from rpy2.robjects.packages import importr
edger = importr("edgeR")

from diff_expr.diffexpr import Limma, remove_batch


def get_attribute_df(filenames, data_folder, attribute):
    all_genes = []
    group_pre_post = pd.Series()

    for file in filenames:
        #sample = re.match('.*OnPrem', file.split('/')[-1]).group()
        sample = '.'.join(file.split('/')[-1].split('.')[:-3])
        if file.split('.')[-1] == 'pickle':
            gc_df = pd.read_pickle(file)
        else:
            out_dir = data_folder+'/'+attribute
            out_downloaded_data = out_dir + '/' + file.split('/')[-1] + '.pickle'
            gc_df = pd.read_csv(file, delimiter='\t', encoding="ISO-8859-1", header=2,
                                index_col=['Name', 'Description'], usecols=['Name', 'Description', sample])

            if not os.path.isdir(out_dir):
                os.makedirs(out_dir)
            with open(out_downloaded_data, 'wb') as f:
                pickle.dump(gc_df, f)
        all_genes.append(gc_df)

        if sample.find('Pre') == -1:  # includes 'Peri' samples, which are all pseudo-post
            group_pre_post[sample] = 'posttreatment'
        else:
            group_pre_post[sample] = 'pretreatment'

    # get rnaseqc2_gene_counts file names***
    all_genes = pd.concat(all_genes, axis=1)

    #all_genes.insert(0,'gene',all_genes.index.get_level_values('Name')+'_'+all_genes.index.get_level_values('Description'))

    return all_genes, group_pre_post


def remove_non_paired_samples(sample_df, separate_pres):
    """Remove samples from dataframe that don't have paired pre/post.
    
    Input dataframe should include sample names as index, participant and pre_post as attributes."""
    if separate_pres:
        rna_method = sample_df.groupby(by=['participant', 'pdb_preservation_method', 'pre_post']).size().unstack(
            fill_value=0)
        rna_ff = rna_method.loc[(slice(None), "FF"), :]
        rna_ffpe = rna_method.loc[(slice(None), "FFPE"), :]
        rna_ff_p = rna_ff[(rna_ff['Post'] > 0) & (rna_ff['Pre'] > 0)].reset_index()['participant'].values
        rna_ffpe_p = rna_ffpe[(rna_ffpe['Post'] > 0) & (rna_ffpe['Pre'] > 0)].reset_index()['participant'].values
        return sample_df[((sample_df['participant'].isin(rna_ff_p)) & (sample_df['pdb_preservation_method'] == 'FF')) |
                         ((sample_df['participant'].isin(rna_ffpe_p)) & (sample_df['pdb_preservation_method'] == 'FFPE'))]
    else:
        pre_post_rna_p = set(sample_df[sample_df['pre_post'] == 'Pre']['participant'].unique()) & \
                         set(sample_df[sample_df['pre_post'] == 'Post']['participant'].unique())
        return sample_df[sample_df['participant'].isin(pre_post_rna_p)].copy()


def df_for_samples(all_counts, all_tpm, keep_samples, metrics_df, delivery_ws, outcomes_fn, samples_fn, participant_table=None):
    """Return gene_counts and gene_tpm for only the specified samples (in specified order)."""
    output_gene_counts = all_counts[keep_samples]
    output_gene_tpm = all_tpm[keep_samples]
    output_groups = df_group_vars(metrics_df, keep_samples, delivery_ws, outcomes_fn, samples_fn, participant_table=participant_table)
    
    return output_gene_counts, output_gene_tpm, output_groups


def df_group_vars(metrics_df, selected_samples, delivery_ws, outcomes_fn, samples_fn, participant_table=None):
    metrics_out = metrics_df.loc[selected_samples]

    # get participant number
    participant_out = pd.Series(metrics_out['participant'].str.extract('GBM.ICB-([0-9]+)', expand=False),
                                selected_samples, name='participant')

    # get (process/delivery) batches
    wm_delivery = dalmatian.WorkspaceManager(delivery_ws)
    sample_sets = wm_delivery.get_sample_sets()  # todo need correct sample sets
    samples = wm_delivery.get_samples()
    batches = sample_sets['samples']

    samples_portal = pd.read_csv(samples_fn, sep='\t', index_col='gp_sample_id', comment='#')
    batch_dict = {}
    for b in list(batches.index):
        try:
            d = {samples_portal.loc[s]['collaborator_sample_id']: b for s in batches.loc[b]}
        except KeyError:
            continue
        else:
            batch_dict.update(d)

    batch_series = pd.Series(batch_dict, name="batch").loc[selected_samples]
    batch_series.replace({b: i for i, b in enumerate(batch_dict.keys())}, inplace=True)

    # get preservation method
    preservation_method = metrics_out['pdb_preservation_method']

    # get pre post status todo should I differentiate Peri?
    pre_post_status = metrics_out['pre_post']

    # get timing from treatment
    timing = metrics_out.apply(lambda x: x['dftx_start'] if x['pre_post'] == 'Pre' else x['dftx_end'], axis=1)
    timing.name = 'df_tx'

    timing_abs = timing.astype(float)
    timing_abs.name = 'df_tx_absolute_pre'
    pre_mask = metrics_out[metrics_out['pre_post'] == 'Pre'].index
    timing_abs.loc[pre_mask] = abs(timing_abs.loc[pre_mask])

    # get outcomes
    outcomes = pd.read_csv(outcomes_fn, sep='\t',
                           index_col='PatientId')
    all_participants = ['GBM.ICB-' + num for num in participant_out.unique()]
    outcomes = outcomes.filter(items=all_participants, axis=0)
    participant_helper = participant_out.to_frame().copy()
    participant_helper['helper'] = ['GBM.ICB-' + num for num in participant_out]
    outcomes_out = outcomes.merge(participant_helper, how='right', right_on='helper', left_index=True)
    desired_outcomes = ['OS_ICB_start_15m', 'OS_ICB_start_18m', 'OS_ICB_start_24m', 'Responder', 'OS_ICB_start_Death-LFU', 'Deceased']
    new_names = ['os_15', 'os_18', 'os_24', 'response', 'os_icb_start', 'deceased']
    outcomes_out = outcomes_out[desired_outcomes].rename(columns={old: new for old, new in zip(desired_outcomes, new_names)})
    # todo add in IDH and MGMT status
    # todo gender??

    # get purity
    purity = metrics_out['wxs_purity']
    purity.name = 'purity'

    # get clinical data from participant_table
    if participant_table is not None:
        participant_out = participant_helper.reset_index().rename(
            columns={'index': 'sample_id'}).merge(
            participant_table, how='left', left_on='helper', right_on='participant_id').set_index('sample_id').join(
            metrics_out['pdb_collection_date_dfd'])
        participant_out['age_at_collection_date'] = participant_out['Age at GBM Dx'].astype(float) + \
                                                    participant_out['pdb_collection_date_dfd'].astype(float).divide(365)
        participant_out.drop(columns=['pdb_collection_date_dfd', 'helper', 'participant_id'], inplace=True)

    return pd.concat([pre_post_status, batch_series, participant_out, preservation_method, purity, timing, timing_abs, outcomes_out], axis=1)


def choose_samples(metrics_df, title, output_dir, blocklist=None, goodlist=None, lt=None, gt=None, best_timing=True,
                   best_qc=None, only_paired=True, best_pre_post=True, separate_pres=False, plot_x=None, plot_y=None):
    """Choose best samples from metrics dataframe according to inputs.
    
    Blocklist removes these samples from contention; goodlist is a convenience (reciprocal of blocklist) to only use these samples. Next priority is lt/gt, given as dictionaries of {attribute: value} (as sample must be less than or greater than this value for this attribute). If best_pre_post is True, best_timing, and best_qc define how to sort the remaining samples (to choose which are indeed the best pre and post sample. Only one Pre/Post will be given per participant, defined by the given attributes. Finally, if only_paired is True, only participants with at least one pre and post sample will be returned (no unmatched pre/post samples).
    
    Returns list of sample names and plot for chosen samples.
    
    Best_qc given as (attribute, bool_ascending)
    """
    metrics_selected = metrics_df.copy()
    if not blocklist:
        blocklist=[]
    if goodlist:
        metrics_selected = metrics_selected.loc[goodlist]
    if not lt:
        lt = {}
    if not gt:
        gt = {}

    metrics_selected.drop(index=set(blocklist) & set(metrics_selected.index), inplace=True)  # remove blocked samples
    
    # remove samples that don't meet metric thresholds
    for att, val in lt.items():
        metrics_selected = metrics_selected[(metrics_selected[att] < val) | (metrics_selected[att].isnull())]
    for att, val in gt.items():
        metrics_selected = metrics_selected[(metrics_selected[att] > val) | (metrics_selected[att].isnull())]
    
    # remove samples that don't have at least one pre and one post
    if only_paired:
        metrics_selected = remove_non_paired_samples(metrics_selected, separate_pres=separate_pres)
    if best_pre_post:  
        # sort remaining samples according to best_timing/qc - only has an effect if best_pre_post is True
        metrics_selected['dftx_end'] = metrics_selected['dftx_end'].abs()
        if best_timing:
            metrics_selected.sort_values(by=['participant', 'pre_post', 'dftx_start', 'dftx_end'], ascending=[True, True, False, True], inplace=True)
        elif best_qc:
            metrics_selected.sort_values(by=['participant', 'pre_post', best_qc[0]], ascending=[True, True, best_qc[1]], inplace=True)

        if separate_pres:
            metrics_selected.drop_duplicates(subset=['participant', 'pre_post', 'pdb_preservation_method'],
                                             keep='first', inplace=True)
        else:
            metrics_selected.drop_duplicates(subset=['participant', 'pre_post'], keep='first', inplace=True)
    
    desired_samples = metrics_selected.index.tolist()
    
    # plotting todo make more complete
    fig = plot_samples(metrics_df, desired_samples, title, output_dir=output_dir, x=plot_x, y=plot_y)
    if output_dir:
        pickle.dump(desired_samples, open(os.path.join(output_dir, f'{title}_chosen_samples.pkl'), 'wb'))

    return desired_samples, fig


def plot_samples(metrics_df, selected_samples, title, x=None, y=None, output_dir=None,
                 col='pre_post', **kwargs):
    metrics_df['selected'] = metrics_df.index.map(lambda s: 'Selected' if s in selected_samples else 'Not Selected')
    metrics_df['participant_int'] = metrics_df['participant'].str.extract('GBM.ICB-([0-9]+)', expand=False)

    x = x if x else 'Median Exon CV'
    y = y if y else 'Genes Detected'

    # define default plotting arguments
    kwargs_display = kwargs.copy()
    if 'facet_col' not in kwargs_display.keys():
        kwargs_display['facet_col'] = 'pre_post'
    if 'color' not in kwargs_display.keys():
        color = list(metrics_df['participant_int'])
    else:
        color = list(metrics_df[kwargs_display['color']])
        kwargs_display.pop('color')
        kwargs.pop('color')
    if 'symbol' not in kwargs_display.keys():
        kwargs_display['symbol'] = 'pdb_preservation_method'
    if 'color_discrete_sequence' not in kwargs_display.keys():
        kwargs_display['color_discrete_sequence'] = px.colors.qualitative.Alphabet
    if 'hover_data' not in kwargs_display.keys():
        kwargs_display['hover_data'] = ['pre_post', 'participant', 'index', 'pdb_preservation_method', 'dftx_start', 'dftx_end', 'selected']
    fig = metrics_df.reset_index().plot.scatter(x, y,
                                                color=color,
                                                backend='plotly',
                                                facet_row='selected',
                                                **kwargs_display,
                                                height=800)

    if output_dir:
        if 'facet_col' not in kwargs.keys():
            kwargs['facet_col'] = 'pre_post'
        fig_static = metrics_df.reset_index().plot.scatter(x, y,
                                                           color=color,
                                                           backend='plotly',
                                                           facet_row='selected',
                                                           **kwargs,
                                                           height=800)
        fig_static.write_image(output_dir + f'{title}_selected_samples.svg')
    return fig


def plot_samples_full_metrics(metrics_df, selected_samples, title, metrics_list=None):
    """fig = px.histogram(sample_metrics, x=metrics_list[3], nbins=30, color='pre_post')
    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.65) """

    if not metrics_list:
        metrics_list = ['Genes Detected', 'Exonic Rate',
                        'Genes Detected', 'Median Exon CV',
                        'Duplicate Rate of Mapped', 'Exonic Reads',
                        'rRNA Rate', "Median 3' bias",
                        'Mapping Rate', 'Low Quality Reads']
    sample_metrics = metrics_df.loc[selected_samples]
    fig = px.scatter_matrix(sample_metrics.reset_index(), dimensions=metrics_list, color='pre_post', title=title,
                            hover_data=['index'])
    fig.update_traces(diagonal_visible=False, showupperhalf=False)
    fig.update_layout(height=800)
    return fig


def get_time_from_tx(participant, sample_dfd, icb_treatments, tx_timing='start'):
    this_tx = icb_treatments[icb_treatments['participant_id'] == participant]
    if this_tx.shape[0] == 0:
        return np.NaN
    # adjusts for one participant (246) that has two ICB treatments todo make more exact
    tx_dfd = this_tx.reset_index().loc[0, tx_timing + '_date_dfd']
    return float(sample_dfd) - float(tx_dfd)


def filter_genes(output_counts, output_tpm, fraction=None, tpm_min=None, count_min=None):
    """Return the indices for the genes to drop and final tpm/counts, based on the given filter values.

    Can possibly extend later to take in a filter function to make more flexible.

    :param fraction: Percent of samples that must meet the minimum requirements for each gene, given as decimal [0,1]
    """
    if not fraction:
        fraction = 0.25
    if not tpm_min:
        tpm_min = 0.1
    if not count_min:
        count_min = 6
    sample_num = output_tpm.shape[1] * fraction
    drop_genes = output_tpm[(output_tpm.gt(tpm_min).sum(axis=1).lt(sample_num)) |
                            (output_counts.gt(count_min).sum(axis=1).lt(sample_num))].index
    dropped_gene_counts = output_counts.drop(index=drop_genes)
    dropped_tpm = output_tpm.drop(index=drop_genes)

    return dropped_gene_counts, dropped_tpm, drop_genes


def compute_correlations(counts, sv_df, batch_group, selected_samples, qc_metrics_df):
    """I want to check the correlation between counts (log counts at least if not log cpm) and a number of batches:
     - All the calculated SVs (stored in self.sv_df of the Limma object)
     - The biological variable (usually pre/post)
     - Any known batch effect (preservation_method) - want to also check batches, but there are more than two categories
     - Days from treatment (absolute value of df_tx for pre only?)
     - QC metrics

     Should return correlation coefficient for every gene for each of these batches (gene x batch matrix)
     Can then turn this into a heatmap
     """
    pertinent_batch = batch_group[['pre_post', 'pdb_preservation_method', 'df_tx', 'os_15', 'os_18', 'os_24', 'response']].copy()
    pertinent_batch['df_tx'] = pertinent_batch.apply(lambda x: abs(float(x['df_tx'])) if x['pre_post'] == 'Pre' else float(x['df_tx']), axis=1)
    if len(np.unique(pertinent_batch['pdb_preservation_method'])) == 1:
        pertinent_batch.drop(columns='pdb_preservation_method', inplace=True)
    full_batch_matrix = pd.concat([pertinent_batch, sv_df, qc_metrics_df], axis=1)
    full_batch_matrix = full_batch_matrix.loc[selected_samples]
    counts = counts[selected_samples]

    return general_correlations(counts, full_batch_matrix)


def general_correlations(x_axis, y_axis):
    corr_matrix = []
    for name, x_arr in x_axis.iterrows():
        corr_arr = []
        for y_name, y_arr in y_axis.transpose().iterrows():
            try:
                nan_vals = np.asarray([pd.isnull(y) for y in y_arr])
            except TypeError:
                nan_vals = np.array([False] * len(y_arr.values))
            x_no_null = x_arr.values[~nan_vals]
            y_no_null = y_arr.values[~nan_vals]
            corr_arr.append(compute_single_corr(x_no_null, y_no_null))
        corr_matrix.append(corr_arr)

    df = pd.DataFrame(corr_matrix, index=x_axis.index, columns=y_axis.columns)
    col_sums = df.abs().sum(axis=0)
    return df.drop(columns=col_sums[col_sums == 0].index)


def compute_single_corr(counts, batch_arr):
    if isinstance(batch_arr[0], float) or (isinstance(batch_arr[0], int) and not set(np.unique(batch_arr)) == {0, 1}):
        return compute_cont_cont_corr(counts, batch_arr)
    else:
        return compute_cat_cont_corr(counts, batch_arr)


def compute_cont_cont_corr(counts, batch_arr):
    corr, p = ss.spearmanr(counts, batch_arr)
    return corr


def compute_cat_cont_corr(counts, batch_arr):
    """Helpful explanation of various correlation measures:
    https://medium.com/@outside2SDs/an-overview-of-correlation-measures-between-categorical-and-continuous-variables-4c7f85610365

    The point biserial correlation measure assumes normality and homoscedascity, although these are already assumed by Limma so there isn't too much of an issue there. And it computs a value between -1 and +1, making it easy to compare to Spearman's.
    """
    if len(np.unique(batch_arr)) != 2:
        # print(f'Input categorical variable must contain only two categories, not {np.unique(batch_arr)}')  # todo change to log
        return 0
    elif not (set(np.unique(batch_arr)) == {0, 1}):
        batch_arr = [v == batch_arr[0] for v in batch_arr]
    corr, p = ss.pointbiserialr(batch_arr, counts)
    return corr


def plot_sv_correlation(limma, output_groups, metrics, selected_samples, title, output_dir, selected_columns=None):
    for c in output_groups:
        if len(output_groups[c].unique()) == 1:
            output_groups.drop(columns=c, inplace=True)

    full_batch_matrix = pd.concat([output_groups, metrics.loc[:, 'Mapping Rate':'Exon CV MAD']], axis=1)

    corr_mat_svs = general_correlations(limma.sv_df.transpose(), full_batch_matrix.loc[selected_samples])
    if selected_columns is not None:
        corr_mat_svs = corr_mat_svs[selected_columns]

    fig_height = 10
    if limma.sv_df.shape[1] < 12:
        fig_height = limma.sv_df.shape[1] / 2 + 4
    fig, ax = plt.subplots(figsize=(20, fig_height))
    im = ax.imshow(corr_mat_svs.values, aspect='auto', interpolation='None', vmin=-1, vmax=1, cmap='seismic')
    plt.xticks(np.arange(0, len(corr_mat_svs.columns)), labels=list(corr_mat_svs.columns), rotation=90)
    cbar = plt.colorbar(im)
    cbar.set_label("Correlation", fontsize=14)
    plt.ylabel("SVs", fontsize=18)
    plt.title(f"{title}: SVs correlated with known covariates", fontsize=20)

    if output_dir:
        plt.savefig(output_dir + f'{title}_sv_corr.svg', bbox_inches='tight')
    return fig, corr_mat_svs


def plot_gene_correlation(limma, output_groups, metrics, selected_samples, title, output_dir, gene_num=None):
    if not gene_num:
        gene_num = limma.results_df[limma.results_df['pval_adj'] < 0.05].shape[0]
    # already sorted by pval_adj by (abs) log2FC
    sig_genes = limma.results_df.iloc[0:gene_num].index
    ensg_names = [re.search('ENSG[0-9]+.[0-9]+', e).group() for e in sig_genes]

    counts_df = pd.DataFrame(limma.adjusted_counts, index=limma.counts_df.index, columns=limma.counts_df.columns)

    for c in output_groups:
        if len(output_groups[c].unique()) == 1:
            output_groups.drop(columns=c, inplace=True)

    # pertinent_batch = batch_group[['pre_post', 'pdb_preservation_method', 'df_tx', 'os_15', 'os_18', 'os_24', 'response']].copy()

    full_batch_matrix = pd.concat([output_groups, limma.sv_df, metrics.loc[:, 'Mapping Rate':'Exon CV MAD']], axis=1).loc[selected_samples]
    corr_only = general_correlations(counts_df.loc[ensg_names, selected_samples], full_batch_matrix)

    fig_height = 10
    if gene_num < 24:
        fig_height = limma.sv_df.shape[1] / 4 + 4
    fig, ax = plt.subplots(figsize=(20, fig_height))
    im = ax.imshow(corr_only.values, aspect='auto', interpolation='None', vmin=-1, vmax=1, cmap='seismic')
    plt.xticks(np.arange(0, len(corr_only.columns)), labels=list(corr_only.columns), rotation=90)
    cbar = plt.colorbar(im)
    cbar.set_label("Correlation", fontsize=14)
    _ = plt.ylabel("Significant Genes", fontsize=18)
    plt.title(f"{title}: Top {gene_num} genes correlated with all covariates", fontsize=20)

    if output_dir:
        plt.savefig(output_dir + f'{title}_gene_corr.svg', bbox_inches='tight')

    return fig, corr_only

# ----- PCA Plots ----------

def plot_pca(limma, covariates, log_before_pca=False, use_defaults=True, output_dir=None, gene_num=None, **kwargs):
    if not gene_num:
        gene_num = limma.results_df[limma.results_df['pval_adj'] < 0.05].shape[0]
    sig_genes = limma.results_df.iloc[0:gene_num].index
    ensg_names = [re.search('ENSG[0-9]+.[0-9]+', e).group() for e in sig_genes]

    covariates = covariates.loc[limma.counts_df.columns].copy()

    # use adjusted counts (log normalized gene counts, adjusted for covariates if present)
    adj_counts_df = pd.DataFrame(limma.adjusted_counts, index=limma.counts_df.index, columns=limma.counts_df.columns)
    x_tr_df, pca = calc_pca(adj_counts_df.loc[ensg_names], log=log_before_pca)

    timing_w_pres = [time + '_' + pres for time, pres in zip(covariates['pre_post'].values,
                                                             covariates['pdb_preservation_method'].values)]
    covariates['timing_w_pres'] = timing_w_pres

    if 'title' not in kwargs.keys():
        title = ''
    else:
        title = kwargs['title']
        kwargs.pop('title')
    if 'x' not in kwargs.keys():
        x_val = 'PC1'
    else:
        x_val = kwargs['x']
        kwargs.pop('x')
    if 'y' not in kwargs.keys():
        y_val = 'PC2'
    else:
        y_val = kwargs['y']
        kwargs.pop('y')

    kwargs_display = kwargs.copy()
    if 'symbol_map' not in kwargs_display.keys() and use_defaults:
        kwargs_display['symbol_map'] = {'Pre_FF': 'circle',
                                'Pre_FFPE': 'circle-open-dot',
                                'Post_FF': 'cross',
                                'Post_FFPE': 'cross-open-dot'}
    if 'symbol' not in kwargs_display.keys() and use_defaults:
        kwargs_display['symbol'] = 'timing_w_pres'
    if 'color' not in kwargs_display.keys() and use_defaults:
        kwargs_display['color'] = 'participant'
    if 'color_discrete_sequence' not in kwargs_display.keys() and use_defaults:
        kwargs_display['color_discrete_sequence'] = px.colors.qualitative.Alphabet
    if 'hover_data' not in kwargs_display.keys() and use_defaults:
        kwargs_display['hover_data'] = ['index', 'pre_post', 'participant', 'pres_type', 'Genes Detected', 'Exonic Rate',
                                        'Median Exon CV', 'Duplicate Rate of Mapped', 'rRNA Rate', "Median 3' bias"]

    fig = scatter_plot_pca(covariates.reset_index(), x_tr_df, pca, x=x_val, y=y_val, title=f"{title} - PCA of top {gene_num} genes", **kwargs_display)
    fig.update_traces(marker={'size': 10})

    if output_dir:
        if 'symbol' not in kwargs.keys() and use_defaults:
            kwargs_display['symbol'] = 'pres_type'
        if 'color' not in kwargs.keys() and use_defaults:
            kwargs_display['color'] = 'pre_post'
        fig_static = scatter_plot_pca(covariates.reset_index(), x_tr_df, pca, x=x_val, y=y_val,
                                   title=f"{title} - PCA of top {gene_num} genes", **kwargs_display)
        fig_static.write_image(output_dir + f'{title}_pca.svg')
        fig.write_html(output_dir + f'{title}_pca.html')
    return fig


def calc_pca(gene_counts, log=False):
    if log:
        gene_counts = np.log10(gene_counts + 0.0001).values
    pca = PCA(n_components=4)
    x_tr = pca.fit_transform(gene_counts.transpose())
    x_tr_df = pd.DataFrame(x_tr, columns=[f'PC{i + 1}' for i in range(x_tr.shape[1])])
    return x_tr_df, pca


def scatter_plot_pca(covariates, x_tr_df, pca, x='PC1', y='PC2', **kwargs):
    assert covariates.shape[0] == x_tr_df.shape[0]
    plot_data = pd.concat([covariates, x_tr_df], axis=1)

    if 'labels' not in kwargs.keys():
        labels = {}
    else:
        labels = kwargs['labels']
        kwargs.pop('labels')

    for axis_val in [x, y]:
        if axis_val[:2] == 'PC':
            pc_num = int(axis_val[2:])
            labels = {axis_val: f'{axis_val} ({pca.explained_variance_ratio_[pc_num - 1] * 100:.1f}%)', **labels}

    fig = px.scatter(plot_data.reset_index(), x, y, labels=labels, **kwargs)
    fig.update_traces(marker={'size': 10, 'opacity': 0.6})
    fig.update_layout(width=900, height=650)
    return fig


def pca_covariate_iterations(limma_sv_runs, metrics, gene_num=200, x='PC1', y='PC2', output_dir=None, **kwargs):
    limma_num = len(limma_sv_runs)
    row_num = math.ceil(limma_num / 4)
    col_num = min(limma_num, 4)
    fig = make_subplots(rows=row_num, cols=col_num, subplot_titles=[f'{i} SVs' for i in range(len(limma_sv_runs))])
    
    if 'title' not in kwargs.keys():
        title = ''
    else:
        title = kwargs['title']
        kwargs.pop('title')

    for i, this_limma in enumerate(limma_sv_runs):
        fig_sub = plot_pca(this_limma, metrics, x=x, y=y, log_before_pca=False,
                           use_defaults=False, gene_num=gene_num,
                           **kwargs)
        this_row = math.ceil((i + 1) / 4)
        this_col = (i % 4) + 1
        fig.add_trace(fig_sub['data'][0], row=this_row, col=this_col)
        fig.update_xaxes(title_text=fig_sub.layout['xaxis_title_text'], row=this_row, col=this_col)
        fig.update_yaxes(title_text=fig_sub.layout['yaxis_title_text'], row=this_row, col=this_col)

    fig.update_layout(height=400 * row_num, title_text=f'PCA on SV iterations - {title}')
    
    if output_dir:
        fig.write_image(os.path.join(output_dir, f'{title}_pca_sv_iterations.svg'))
    return fig

# -------------------------------

def run_full_analysis(metrics, all_gene_counts, all_gene_tpm, outcomes_fn, samples_fn, title, output_dir=None,
                      delivery_ws='terra-broad-cancer-prod/Getz_Wu_IBM_GBM_RNAseq', gene_qc_fraction=None,
                      gene_tpm_min=None, gene_count_min=None, limma_run_sva=True, limma_nsv=None, limma_att=None,
                      limma_order=None, limma_gene_name_s=None, limma_contrast_def=None, limma_coef=None,
                      limma_model=None, voom_quality_weights=False, corr_gene_num=None, pca_gene_num=None,
                      correlation_id=None, participant_table=None,
                      **samples_att):
    if limma_order is None:
        limma_order = ['Pre', 'Post']
    if limma_att is None:
        limma_att = ['participant', 'pre_post']  # condition should be last

    selected_samples, samples_fig = choose_samples(metrics, title, output_dir, **samples_att)
    output_gene_counts, output_gene_tpm, output_groups = df_for_samples(all_gene_counts, all_gene_tpm,
                                                                        selected_samples, metrics,
                                                                        delivery_ws=delivery_ws,
                                                                        outcomes_fn=outcomes_fn,
                                                                        samples_fn=samples_fn,
                                                                        participant_table=participant_table)
    final_gene_counts, final_tpm, dropped_genes = filter_genes(output_gene_counts, output_gene_tpm,
                                                               fraction=gene_qc_fraction, tpm_min=gene_tpm_min,
                                                               count_min=gene_count_min)
    limma = Limma(output_gene_counts, output_groups[limma_att], limma_order, dropped_genes,
                  gene_name_s=limma_gene_name_s)
    limma.run(run_sva=limma_run_sva, nsv=limma_nsv, contrast_def=limma_contrast_def,
              coef=limma_coef, model_def=limma_model, voom_quality_weights=voom_quality_weights,
              correlation_id=correlation_id)

    # call plotting utilities
    correlation_output = output_groups.drop(columns=['deceased', 'os_icb_start', 'batch', 'participant'])
    gene_corr_fig, _ = plot_gene_correlation(limma, correlation_output, metrics, selected_samples, title, output_dir,
                                             gene_num=corr_gene_num)
    if limma_run_sva:
        sva_corr_fig, _ = plot_sv_correlation(limma, correlation_output, metrics, selected_samples, title, output_dir)
    else:
        sva_corr_fig = None
    pca_fig = plot_pca(limma, metrics, selected_samples, title, output_dir, gene_num=pca_gene_num)

    # return plot named tuple and named tuple for other data
    fig_dict = {'selected_samples': samples_fig,
                'gene_correlation': gene_corr_fig,
                'sva_correlation': sva_corr_fig,
                'gene_pca': pca_fig,
                'voom_mean_variance': Image(data=limma.voom_plot_bytes, format='jpeg', embed=True)}
    results_dict = {'selected_samples': selected_samples,
                    'gene_counts': final_gene_counts,
                    'output_groups': output_groups,
                    'limma': limma}

    return results_dict, fig_dict


def run_full_rna_iteration(output_gene_counts, design, dropped_genes, title, max_nsv=None, outdir=None, run_spec='',
                           coef='os_icb_start', biological_model='os_icb_start', model_def='os_icb_start'):
    # calculate max_nsv programmatically
    if max_nsv is not None:
        nsv_dict = {'nsv': max_nsv}
    else:
        nsv_dict = {}
    full_sv_limma = Limma(output_gene_counts, design, None, dropped_genes)
    full_sv_limma.run(run_sva=True, voom_quality_weights=True,
                      coef=coef, biological_model=biological_model,
                      model_def=model_def, **nsv_dict)

    sv_num = full_sv_limma.sv_df.shape[1]
    store_runs = []
    store_gsea = []
    full_gsea_df = []
    full_de_genes_df = []
    for i in range(sv_num + 1):
        if i == 0:
            this_limma = Limma(output_gene_counts, design, None, dropped_genes)
            this_limma.run(run_sva=False, voom_quality_weights=True,
                           coef=coef, biological_model=biological_model,
                           model_def=model_def)
        else:
            design_mod = design.join(full_sv_limma.sv_df.loc[:, :f'SV{i}'])
            this_limma = Limma(output_gene_counts, design_mod, None, dropped_genes)
            this_limma.run(run_sva=False, voom_quality_weights=True,
                           coef=coef, biological_model=biological_model,
                           model_def=' + '.join(list(design_mod)))
            this_limma.adjusted_counts = remove_batch(this_limma.adjusted_counts,
                                                      covariates=this_limma.mod.iloc[:, 2:],
                                                      design=this_limma.mod.iloc[:, :2])

        this_limma.qq_plot(title=f'{title} - Num SVs: {i}')
        store_runs.append(this_limma)
        gsea_res = this_limma.run_gsea()
        store_gsea.append(gsea_res)
        gsea_res.dotplot()
        results = gsea_res.get_results()
        results['iteration'] = i
        full_gsea_df.append(results)

        results_genes = this_limma.results_df
        results_genes['iteration'] = i
        full_de_genes_df.append(results_genes)

    full_gsea_df = pd.concat(full_gsea_df)
    full_de_genes_df = pd.concat(full_de_genes_df)

    full_de_genes_df.reset_index(inplace=True)
    full_de_genes_df['ensg_gene'] = full_de_genes_df['gene_id'].apply(lambda x: x.split("', '")[0][2:])
    full_de_genes_df['hugo_symbol'] = full_de_genes_df['gene_id'].apply(lambda x: x.split("', '")[1][:-2])

    if outdir is not None:
        save_full_rna_iterations(outdir, run_spec, store_runs, store_gsea,
                                 full_gsea_df, full_de_genes_df, full_sv_limma)

    return store_runs, store_gsea, full_gsea_df, full_de_genes_df, full_sv_limma


def save_full_rna_iterations(outdir, run_spec, store_runs, store_gsea, full_gsea_df, full_de_genes_df, sv_limma_or_samples, kfold=False):
    pickle.dump(store_runs, open(os.path.join(outdir, run_spec + 'all_limma_runs.pkl'), 'wb'))
    pickle.dump(store_gsea, open(os.path.join(outdir, run_spec + 'all_gsea_runs.pkl'), 'wb'))
    pickle.dump(full_gsea_df, open(os.path.join(outdir, run_spec + 'gsea_df.pkl'), 'wb'))
    pickle.dump(full_de_genes_df, open(os.path.join(outdir, run_spec + 'all_de_genes.pkl'), 'wb'))
    
    fn = 'samples.pkl' if kfold else 'full_sv_limma.pkl'
    pickle.dump(sv_limma_or_samples, open(os.path.join(outdir, run_spec + fn), 'wb'))


def load_full_rna_iterations(outdir, run_spec, kfold=False):
    store_runs = pd.read_pickle(os.path.join(outdir, run_spec + 'all_limma_runs.pkl'))
    store_gsea = pd.read_pickle(os.path.join(outdir, run_spec + 'all_gsea_runs.pkl'))
    full_gsea_df = pd.read_pickle(os.path.join(outdir, run_spec + 'gsea_df.pkl'))
    full_de_genes_df = pd.read_pickle(os.path.join(outdir, run_spec + 'all_de_genes.pkl'))
    
    fn = 'samples.pkl' if kfold else 'full_sv_limma.pkl'
    sv_limma_or_samples = pd.read_pickle(os.path.join(outdir, run_spec + fn))

    return store_runs, store_gsea, full_gsea_df, full_de_genes_df, sv_limma_or_samples


def plot_gsea_incr_sv(full_gsea_df, qval_cutoff=0.15, only_plot_sig=None, output_dir=None, run_spec=''):
    full_gsea_df['logFDR'] = full_gsea_df['FDR q-val'].apply(lambda x: np.abs(np.log(x + 0.0000001)))
    full_gsea_df['dir'] = full_gsea_df['NES'].gt(0).replace({True: 'Postive NES', False: 'Negative NES'})

    only_sig_gsea_df = full_gsea_df[full_gsea_df['FDR q-val'] < qval_cutoff]
    flipping_sig_gsea = only_sig_gsea_df.groupby(['Term', 'dir']).size().unstack().dropna()
    flipping_gsea = full_gsea_df.groupby(['Term', 'dir']).size().unstack().dropna()

    never_sig = set(full_gsea_df['Term']) - set(only_sig_gsea_df['Term'])
    no_flips_sig = set(only_sig_gsea_df['Term']) - set(flipping_gsea.index)
    flips_sig_both = set(flipping_sig_gsea.index)
    flips_sig_one = set(only_sig_gsea_df['Term']) - flips_sig_both - no_flips_sig
    type_dict = dict(**dict.fromkeys(never_sig, 'Never Significant'),
                     **dict.fromkeys(no_flips_sig, 'Significant, No Flips'),
                     **dict.fromkeys(flips_sig_both, 'Significant both Directions'),
                     **dict.fromkeys(flips_sig_one, 'Significant one Direction, but flips'))

    if only_plot_sig is not None:
        only_these_terms = only_sig_gsea_df[only_sig_gsea_df['iteration'] == only_plot_sig]['Term'].values
    else:
        only_these_terms = full_gsea_df['Term'].unique()

    full_gsea_df['Type'] = full_gsea_df['Term'].map(type_dict)
    fig = px.line(full_gsea_df[full_gsea_df['Term'].isin(only_these_terms)], x='iteration', y='NES', color='Term', facet_col='Type', title=f'Increasing SVs - {run_spec}',
                  height=400, width=1200, hover_data=['FDR q-val'])
    # add stars for significance - fig.add_trace(go.Scatter(
    
    if output_dir:
        run_spec = f'{run_spec}_' if run_spec != '' else run_spec
        qval_spec = f'_qval_{qval_cutoff}'
        fig.write_image(os.path.join(output_dir, f'{run_spec}increasing_svs{qval_spec}.svg'))
        fig.write_html(os.path.join(output_dir, f'{run_spec}increasing_svs{qval_spec}.html'))
    return fig


def nes_dir_table(full_gsea_df, filter_on=None, filter_val=0.15, output_dir=None, run_spec=''):
    full_gsea_df['NES_pos'] = full_gsea_df['NES'].gt(0)
    num_pos = full_gsea_df.groupby('Term')['NES_pos'].sum()
    min_qval = full_gsea_df[full_gsea_df['iteration'] > 1].groupby('Term')['FDR q-val'].min()
    min_qval.name = 'Minimum q-val'
    last_qval = full_gsea_df[full_gsea_df['iteration'] == full_gsea_df['iteration'].max()][
        ['Term', 'FDR q-val']].set_index('Term')
    last_qval.rename(columns={'FDR q-val': 'FDR q-val All SVs'}, inplace=True)

    tmp = pd.concat([num_pos, min_qval, last_qval], axis=1)
    if filter_on == 'min':
        mask = tmp['Minimum q-val'] <= filter_val
    elif filter_on == 'final':
        mask = tmp['FDR q-val All SVs'] <= filter_val
    else:
        mask = pd.Series(True, index=tmp.index)
    styled_df = tmp[mask].sort_values(['NES_pos', 'FDR q-val All SVs'],
                                      ascending=[True, False]).style.background_gradient(subset=['NES_pos'],
                                                                            cmap=plt.cm.get_cmap('coolwarm'), vmin=0,
                                                                            vmax=tmp['NES_pos'].max())
    if output_dir:
        filter_on_spec = f'_filter_on_{filter_on}' if filter_on is not None else ''
        filter_val_spec = f'_lt_{filter_val}'
        run_spec = f'{run_spec}_' if run_spec != '' else run_spec
        styled_df.to_html(os.path.join(output_dir, f'{run_spec}gsea_sig_results{filter_on_spec}{filter_val_spec}.html'))
                           
    return styled_df


def run_kfold_rna_iteration(output_gene_counts, design, dropped_genes, title, k_folds=5, max_nsv=None, outdir=None, run_spec='',
                           coef='os_icb_start', biological_model='os_icb_start', model_def='os_icb_start', gsea_dict=None):
    # calculate max_nsv programmatically
    if max_nsv is not None:
        nsv_dict = {'nsv': max_nsv}
    else:
        nsv_dict = {}
        
    if gsea_dict is None:
        gsea_dict = {}
        
    from sklearn.model_selection import KFold
    # split dataset into k folds
    assert k_folds <= design.shape[0]
    kf = KFold(n_splits=k_folds, shuffle=True)
    store_samples = []
    store_sv_num = []
    store_runs = []
    store_gsea = []
    full_gsea_df = []
    full_de_genes_df = []
    for i, (train_index, _) in enumerate(kf.split(design)):
        fold_design = design.iloc[train_index]
        fold_gene_counts = output_gene_counts.loc[:, fold_design.index]
        store_samples.append(fold_design.index)
        
        fold_sv_limma = Limma(fold_gene_counts, fold_design, None, dropped_genes)
        fold_sv_limma.run(run_sva=True, voom_quality_weights=True,
                          coef=coef, biological_model=biological_model,
                          model_def=model_def, **nsv_dict)

        store_sv_num.append(fold_sv_limma.sv_df.shape[1])

        fold_sv_limma.qq_plot(title=f'{title} - Fold: {i}')
        store_runs.append(fold_sv_limma)
        gsea_res = fold_sv_limma.run_gsea(**gsea_dict)
        store_gsea.append(gsea_res)
        gsea_res.dotplot()
        results = gsea_res.get_results()
        results['iteration'] = i
        full_gsea_df.append(results)

        results_genes = fold_sv_limma.results_df
        results_genes['iteration'] = i
        full_de_genes_df.append(results_genes)

    full_gsea_df = pd.concat(full_gsea_df)
    full_de_genes_df = pd.concat(full_de_genes_df)

    full_de_genes_df.reset_index(inplace=True)
    full_de_genes_df['ensg_gene'] = full_de_genes_df['gene_id'].apply(lambda x: x.split("', '")[0][2:])
    full_de_genes_df['hugo_symbol'] = full_de_genes_df['gene_id'].apply(lambda x: x.split("', '")[1][:-2])

    if outdir is not None:
        save_full_rna_iterations(outdir, run_spec, store_runs, store_gsea,
                                 full_gsea_df, full_de_genes_df, store_samples, kfold=True)

    return store_runs, store_gsea, full_gsea_df, full_de_genes_df, store_samples, store_sv_num


def run_bootstrap_rna(output_gene_counts, design, dropped_genes, title, iterations=100, percentile=0.9, leave_out=0.1, max_nsv=None, outdir=None, run_spec='',
                           coef='os_icb_start', biological_model='os_icb_start', model_def='os_icb_start', gsea_dict=None):
    # calculate max_nsv programmatically
    if max_nsv is not None:
        nsv_dict = {'nsv': max_nsv}
    else:
        nsv_dict = {}
        
    if gsea_dict is None:
        gsea_dict = {}
    
    full_sv_limma = Limma(output_gene_counts, design, None, dropped_genes)
    full_sv_limma.run(run_sva=True, voom_quality_weights=True,
                      coef=coef, biological_model=biological_model,
                      model_def=model_def, **nsv_dict)
    
    sample_num = int(np.floor(len(design.index)*(1-leave_out)))
    full_gsea_df = []
    store_samples = []
    for i in range(iterations):
        print(f'Iteration: {i}')
        bootstrap_sample = np.random.choice(design.index, size=sample_num, replace=False)
        store_samples.append(bootstrap_sample)
        design_mod = design.join(full_sv_limma.sv_df).loc[bootstrap_sample]
        this_limma = Limma(output_gene_counts.loc[:, bootstrap_sample], design_mod, None, dropped_genes)
        this_limma.run(run_sva=False, voom_quality_weights=True,
                       coef=coef, biological_model=biological_model,
                       model_def=' + '.join(list(design_mod)))
        this_gsea = this_limma.run_gsea(**gsea_dict).get_results()
        this_gsea['iter'] = i
        full_gsea_df.append(this_gsea)
        
    full_gsea_df = pd.concat(full_gsea_df)
    nes_df = full_gsea_df.set_index(['iter', 'Term'])['NES'].unstack(1)
    ci_low = 50 * (1 - percentile)
    ci_high = 50 * (1 + percentile)
    lower_ci = np.percentile(nes_df.values, ci_low, axis=0)
    upper_ci = np.percentile(nes_df.values, ci_high, axis=0)

    summary_df = pd.DataFrame([lower_ci, upper_ci, nes_df.mean(), nes_df.std()], columns=nes_df.columns, index=[f'lower_CI_{ci_low:.1f}', f'upper_CI_{ci_high:.1f}', 'mean', 'std']).T
    summary_df['same_dir'] = (summary_df[f'lower_CI_{ci_low:.1f}'] * summary_df[f'upper_CI_{ci_high:.1f}']) > 0
    
    full_gsea = full_sv_limma.run_gsea().get_results().set_index('Term')
    summary_df = summary_df.join(full_gsea[['NES', 'NOM p-val', 'FDR q-val']])
    
    # save full_gsea_df
    # save summary
    if outdir is not None:
        pickle.dump(full_gsea_df, open(os.path.join(outdir, run_spec + '_all_gsea_iterations.pkl'), 'wb'))
        pickle.dump(summary_df, open(os.path.join(outdir, run_spec + '_gsea_summary.pkl'), 'wb'))
        summary_df.to_csv(os.path.join(outdir, run_spec + '_gsea_summary.txt'), sep='\t')
        pickle.dump(store_samples, open(os.path.join(outdir, run_spec + '_samples.pkl'), 'wb'))
        pickle.dump(full_sv_limma, open(os.path.join(outdir, run_spec + '_full_sv_limma.pkl'), 'wb'))

    return full_gsea_df, summary_df, store_samples, full_sv_limma


def load_bootstrap_rna(outdir, run_spec):
    full_gsea_df = pd.read_pickle(os.path.join(outdir, run_spec + '_all_gsea_iterations.pkl'))
    summary_df = pd.read_pickle(os.path.join(outdir, run_spec + '_gsea_summary.pkl'))
    store_samples = pd.read_pickle(os.path.join(outdir, run_spec + '_samples.pkl'))
    full_sv_limma = pd.read_pickle(os.path.join(outdir, run_spec + '_full_sv_limma.pkl'))

    return full_gsea_df, summary_df, store_samples, full_sv_limma


def run_leave_k_out_rna(output_gene_counts, design, dropped_genes, title, leave_out=1, max_nsv=None, outdir=None, run_spec='',
                           coef='os_icb_start', biological_model='os_icb_start', model_def='os_icb_start', gsea_dict=None, **kwargs):
    # calculate max_nsv programmatically
    if max_nsv is not None:
        nsv_dict = {'nsv': max_nsv}
    else:
        nsv_dict = {}
        
    if gsea_dict is None:
        gsea_dict = {}
    
    full_sv_limma = Limma(output_gene_counts, design, None, dropped_genes)
    full_sv_limma.run(run_sva=True, voom_quality_weights=True,
                      coef=coef, biological_model=biological_model,
                      model_def=model_def, **nsv_dict)
    
    all_samples = design.index
    combos = itertools.combinations(all_samples, len(all_samples) - leave_out)
    
    full_gsea_df = []
    store_samples = []
    for i, bs in enumerate(combos):
        bootstrap_sample = list(bs)
        print(f'Iteration: {i}')
        store_samples.append(bootstrap_sample)
        design_mod = design.join(full_sv_limma.sv_df).loc[bootstrap_sample]
        this_limma = Limma(output_gene_counts.loc[:, bootstrap_sample], design_mod, None, dropped_genes)
        this_limma.run(run_sva=False, voom_quality_weights=True,
                       coef=coef, biological_model=biological_model,
                       model_def=' + '.join(list(design_mod)))
        this_gsea = this_limma.run_gsea(**gsea_dict).get_results()
        this_gsea['iter'] = i
        full_gsea_df.append(this_gsea)
        
    full_gsea_df = pd.concat(full_gsea_df)
    
    if outdir is not None:
        pickle.dump(full_gsea_df, open(os.path.join(outdir, run_spec + '_all_gsea_iterations.pkl'), 'wb'))
        pickle.dump(store_samples, open(os.path.join(outdir, run_spec + '_samples.pkl'), 'wb'))
        pickle.dump(full_sv_limma, open(os.path.join(outdir, run_spec + '_full_sv_limma.pkl'), 'wb'))
    
    nes_df = full_gsea_df.set_index(['iter', 'Term'])['NES'].unstack(1)
    pval_df = full_gsea_df.set_index(['iter', 'Term'])['NOM p-val'].unstack(1)
    store_samples_missing = [list(set(all_samples) - set(l)) for l in store_samples]
    max_p_values = pval_df.astype(float).max()
    tmp_max = pval_df.astype(float).idxmax()
    max_missing_samples = tmp_max.apply(lambda x: store_samples_missing[x])
    tmp_max_idx_all = pd.Series({k: pval_df[pval_df[k] == v].index.tolist() for k, v in max_p_values.to_dict().items()})
    max_missing_samples_all = tmp_max_idx_all.apply(lambda x: [store_samples_missing[v] for v in x])
    max_nes_values_all = pd.Series({k: [nes_df.loc[v, k] for v in l] for k, l in tmp_max_idx_all.to_dict().items()})
    min_p_values = pval_df.astype(float).min()
    tmp_min_idx_all = pd.Series({k: pval_df[pval_df[k] == v].index.tolist() for k, v in min_p_values.to_dict().items()})
    min_nes_values_all = pd.Series({k: [nes_df.loc[v, k] for v in l] for k, l in tmp_min_idx_all.to_dict().items()})

    summary_df = pd.DataFrame([min_p_values, max_p_values, max_missing_samples, max_missing_samples_all, max_nes_values_all, min_nes_values_all], columns=pval_df.columns, index=['min_pval', 'max_pval', 'missing_samples_max_first', 'missing_samples_max_all', 'max_nes_values_all', 'min_nes_values_all']).T
    
    full_gsea = full_sv_limma.run_gsea(**gsea_dict).get_results().set_index('Term')
    summary_df = summary_df.join(full_gsea[['NES', 'NOM p-val', 'FDR q-val']])
    
    # save data
    if outdir is not None:
        pickle.dump(summary_df, open(os.path.join(outdir, run_spec + '_gsea_summary.pkl'), 'wb'))
        summary_df.to_csv(os.path.join(outdir, run_spec + '_gsea_summary.txt'), sep='\t')

    return full_gsea_df, summary_df, store_samples, full_sv_limma


def plot_comparison_results(kfold_df_x, kfold_df_y, compare_df, display_values=None):
    if display_values is None:
        display_values = compare_df['Term'].values
    
    # get filtered x-axis data (NES and CI values)
    x_nes = kfold_df_x.groupby(['Term'])['NES'].mean().loc[display_values]
    x_ci_max = kfold_df_x.groupby(['Term'])['NES'].max().loc[display_values]
    x_ci_min = kfold_df_x.groupby(['Term'])['NES'].min().loc[display_values]
    
    # get filtered y-axis data (NES and CI values)
    y_nes = kfold_df_y.groupby(['Term'])['NES'].mean().loc[display_values]
    y_ci_max = kfold_df_y.groupby(['Term'])['NES'].max().loc[display_values]
    y_ci_min = kfold_df_y.groupby(['Term'])['NES'].min().loc[display_values]
    
    # get filtered significance values
    sig_values = compare_df.set_index('Term').loc[display_values, ['FDR q-val']]
    sig_values['logFDR'] = sig_values['FDR q-val'].apply(lambda x: np.abs(np.log(x + 0.0000001)))
    sig_values['sig_level'] = np.log10(1 / sig_values['FDR q-val'].astype(float)).replace({np.inf: 3.5}).apply(math.floor)
    # significance color
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x_nes, y=y_nes,
        mode='markers',
        marker=dict(size=sig_values['sig_level'],
                   # color=
                   ),
        error_y=dict(
            type='data',
            symmetric=False,
            array=y_ci_max - y_nes,
            arrayminus=y_nes - y_ci_min
        ),
        error_x=dict(
            type='data',
            symmetric=False,
            array=x_ci_max - x_nes,
            arrayminus=x_nes - x_ci_min
        )))
    
    return fig