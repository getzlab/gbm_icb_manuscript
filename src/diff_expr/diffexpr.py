"""diffexpr.py: Python wrapper classes for R differential expression packages"""

__author__ = "Francois Aguet"
__copyright__ = "Copyright 2020, The Broad Institute"
__license__ = "BSD3"

import pandas as pd
import numpy as np
import subprocess
import copy
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import qtl.plot
import gseapy
import blitzgsea as blitz


has_rpy2 = False
e = subprocess.call('which R', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
try:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.lib import grdevices
    pandas2ri.activate()
    to_dataframe = ro.r('function(x) data.frame(x)')

    limma = importr("limma")
    edger = importr("edgeR")
    deseq2 = importr("DESeq2")
    biomart = importr("biomaRt")
    fgsea = importr("fgsea")

    if e == 0:
        has_rpy2 = True
except:
    pass
if not has_rpy2:
    print("Warning: 'rfunc' cannot be imported. R and the 'rpy2' Python package are needed.")


class Limma(object):
    def __init__(self, full_counts_df, design, order, dropped_genes, gene_name_s=None):
        """
        :param full_counts_df: counts matrix, all_genes x samples
        :param design:  class label(s) for each sample; if running default model, parameter of interest must be in last column
        :param order: test order for the variable of interest; the first element is the reference
        :param dropped_genes: list of gene names to filter out after calculating DGE Norm Factors
        :param gene_name_s: alternate gene names to join to the results dataframe (default None)
        """
        self.method = 'limma-voom'
        if order is not None:
            self.design_df = design.astype('category')  # need to adjust this todo make automatic but flexible
            if isinstance(self.design_df, pd.Series):
                self.design_df = self.design_df.to_frame()  # intercept will be automatically added by model definition

            var_name = self.design_df.columns[-1]
            assert self.design_df[var_name].isin(order).all()
            self.design_df[var_name] = self.design_df[var_name].cat.reorder_categories(order)
        else:
            self.design_df = design

        self.gene_name_s = gene_name_s
        dge_init = edger.DGEList(counts=full_counts_df)
        dge_init = edger.calcNormFactors(dge_init)  # TMM

        self.counts_df = full_counts_df.loc[full_counts_df.index.difference(dropped_genes, sort=False), design.index]
        self.dge = edger.DGEList(counts=self.counts_df)
        self.dge[1] = dge_init[1]

        self.v = None
        self.mod = None
        self.fit = None
        self.contrasts = None
        self.voom_plot_bytes = None
        self.sv_df = pd.DataFrame()
        self.adjusted_counts = edger.cpm(self.dge, log=True)  # adjusted for gene length, normalization and (if SVA is True) covariates

    def run(self, run_sva=False, nsv=None, coef=None, contrast_def=None, model_def=None,
            voom_quality_weights=False, biological_model=None, correlation_id=None):
        """Run voom normalization, fit limma's linear model to data, and calculate contrast based on inputs.

        :param run_sva: Boolean to run SVA batch correction (default: False)
        :param nsv: Number of surrogate variables to use in SVA (default: estimated from data)
        :param coef: Name of coefficient to run differential expression over (default: last coefficient in model)
        :param contrast_def: Define contrasts by combination of model columns
        :param model_def: Define model from design df (default: last column in design)
        :param voom_quality_weights: Boolean to run Voom with quality weights,
            including sample-specific quality weights to increase robostness (default: False)
        :param biological_model: Define biological factor(s) to be protected from SVA;
            must be defined if there is more than 1 factor or if not in final column of design
            (default: last column in design)
        :param correlation_id: Column in design_df to use to compute duplicateCorrelation (default: None)
        :return: None
        """

        # define model matrix (can this be simplified using pd.get_dummies?)
        ro.r.assign('design.df', self.design_df)
        if biological_model is None:  # must be defined if there is more than one factor or if it is not in last column
            biological_model = self.design_df.columns[-1]
        mod = ro.r(f"mod <- model.matrix(as.formula('~{model_def if model_def else biological_model}'), design.df)")
        mod_sv = ro.r(f"mod_sv <- model.matrix(as.formula('~ {biological_model}'), design.df)")

        self.mod = pd.DataFrame(mod, index=ro.r("rownames(mod)"), columns=ro.r('colnames(mod)'))
        with grdevices.render_to_bytesio(grdevices.jpeg, width=512, height=512, res=150) as img:
            if voom_quality_weights and not run_sva:
                v = limma.voomWithQualityWeights(self.dge, design=self.mod, plot=True)  # voom-transformed counts
            else:
                v = limma.voom(self.dge, design=self.mod, plot=True)
        vd = dict(zip(v.names, list(v)))

        if run_sva:
            smartsva = importr("SmartSVA")
            print('Running SmartSVA')
            if nsv is None:  # estimate number of SVs
                isva = importr("isva")
                ro.r.assign('v', v)
                Yr = ro.r(f"t(resid(lm(as.formula('t(v$E) ~ {biological_model}'), data=design.df)))")
                res = isva.EstDimRMT(Yr, False)
                res = dict(zip(res.names, list(res)))
                nsv = res['dim'][0] + 1
            print(f'  * SVs: {nsv}')
            # smartsva: mod0 default is intercept -> mod_sv must have an intercept
            sv_obj = smartsva.smartsva_cpp(vd['E'], mod_sv,
                                           mod0=ro.NULL, n_sv=nsv, alpha=1, B=200, VERBOSE=True)
            sv_obj = dict(zip(sv_obj.names, list(sv_obj)))
            self.sv_df = pd.DataFrame(sv_obj['sv'], index=self.design_df.index,
                                      columns=[f'SV{i}' for i in range(1, sv_obj['sv'].shape[1]+1)])
            # update model with SVs
            self.mod = pd.concat([self.mod, self.sv_df], axis=1)
            with grdevices.render_to_bytesio(grdevices.jpeg, width=512, height=512, res=150) as img:
                if voom_quality_weights:
                    v = limma.voomWithQualityWeights(self.dge, design=self.mod, plot=True)
                else:
                    v = limma.voom(self.dge, design=self.mod, plot=True)

            # mod does not include sv's, only pre_post and perhaps blocking variable
            self.adjusted_counts = self.remove_batch(covariates=self.sv_df, design=mod)

        # fixes issue when using complex multieffect model
        self.mod.columns = [v.replace(':', '_') for v in self.mod.columns]

        self.voom_plot_bytes = img.getvalue()
        self.v = v
        if correlation_id:
            ro.r.assign('mod', self.mod)
            ro.r.assign('v', v)
            id_values = self.design_df[correlation_id].values.astype(int)
            ro.r.assign('block_id', id_values)

            cor = ro.r(f"cor <- duplicateCorrelation(v, mod, block=block_id)")
            print(f'Correlation: {cor[0][0]}')

            self.fit = limma.lmFit(v, design=self.mod, block=id_values, correlation=cor[0][0])
        else:
            self.fit = limma.lmFit(v, design=self.mod)

        if contrast_def:
            ro.r.assign('mod', self.mod)
            contrasts = ro.r(
                f"contrasts <- makeContrasts({contrast_def}, levels=mod)")
            ro.r.assign('fit', self.fit)
            self.fit = ro.r(f"fit <- contrasts.fit(fit, contrasts)")
            self.contrasts = contrasts

        self.fit = limma.eBayes(self.fit)

        _ = self.get_results(coef=coef, save_results=True)

    def get_results(self, coef=None, save_results=True):
        """Calculate contrast for given coefficient

        :param coef: Name of column or contrast to use for DE analysis (if not given, uses last column in model)
        :param save_results: Boolean on whether to save results to self.results_df (default: True)
        :return: pandas.DataFrame giving DE results
        """
        if not self.fit:
            print('Model has not been fit yet. Before accessing results, must run() on Limma object.')
        with localconverter(ro.default_converter + pandas2ri.converter):
            if not coef:
                coef = self.mod.shape[1] - self.sv_df.shape[1]
                print(f'Using Contrast Coefficient from default model (w/ intercept, group variable defined last)')
            print(f'Coefficient: {coef}')
            results_df = ro.conversion.rpy2py(limma.topTable(self.fit, n=np.Inf, coef=coef, sort_by='none'))

        # limma_df.rename(columns={i:i.replace('.','_') for i in limma_df.columns}, inplace=True)
        results_df.rename(columns={'logFC':'log2FC', 'AveExpr':'log2_mean_expr', 't':'tstat',
                                   'P.Value':'pval', 'adj.P.Val':'pval_adj'}, inplace=True)
        results_df.index.name = 'gene_id'

        if self.gene_name_s is not None:
            results_df = self.gene_name_s.rename('gene_name').to_frame().join(results_df, how='right')

        # sort first by pval_adj and then by (absolute_val of) log2FC
        results_df = results_df.assign(A=results_df['log2FC'].abs()).sort_values(['pval_adj', 'A'], ascending=[True, False]).drop('A', axis=1)

        if save_results:
            self.results_df = results_df
            print('Saving results')
        print(f"DE genes at 0.05 FDR: {sum(results_df['pval_adj']<=0.05)}")
        return results_df

    def remove_batch(self, covariates, design):
        return limma.removeBatchEffect(self.adjusted_counts, covariates=covariates, design=design)

    def qq_plot(self, title=''):
        qtl.plot.qqplot(self.results_df['pval'], fontsize=12, title=title)

    def ma_plot(self):
        ma_plot(self.results_df['log2_mean_expr'], self.results_df['log2FC'], gridsize=100)

    def volcano_plot(self):
        volcano_plot(self.results_df['log2FC'],
                     self.results_df['pval_adj'],
                     ylabel='-log$\mathregular{_{10}}$(adj. p-value)')

    def run_gsea(self, gene_set='MSigDB_Hallmark_2020', order_by='sign_pval', **kwargs):
        """Run GSEA prerank on differential expression results

        :param gene_set: Gene set collection to pass to GSEA (default: MSigDB_Hallmark_2020)
        :param order_by: Sort genes by one of sign_pval (default) or weight_pval_logfc
        :param kwargs: other named arguments to pass to gseapy's prerank module
        :return:
        """
        df = self.results_df.reset_index().copy()
        # df['ensembl_gene_id'] = df['gene_id'].apply(lambda x: x.split("'")[1].split(".")[0]) # not used

        df['genes_hugo'] = df['gene_id'].apply(lambda x: x.split("'")[3])
        if order_by == 'sign_pval':
            df['enrichment_score'] = (2 * df['log2FC'].gt(0) - 1) * -np.log10(df['pval'])  # signed p-value
        elif order_by == 'weight_pval_logfc':
            df['enrichment_score'] = df['log2FC'] * 100 * -np.log10(df['pval'])  # p-value weighted by logFC
        else:
            raise(ValueError(f'order_by input must be one of [sign_pval, weight_pval_logfc], not {order_by}'))

        return GSEA(df[['genes_hugo', 'enrichment_score']], gene_set, **kwargs)

    ####################################################
    # Should move to the GSEA object and fully implement
    # def run_camera(self, gene_sets_file='https://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata', inter_gene_cor=0.01):
    #     """CAUTION: Not fully tested yet."""
    #     ensembl = biomart.useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    #     results = biomart.getBM(attributes=['ensembl_gene_id', 'entrezgene_id', 'gene_biotype'],
    #                             filters='ensembl_gene_id',
    #                             values=self.counts_df.index.get_level_values(0),
    #                             mart=ensembl)

    #     readRDS = ro.r['readRDS']
    #     url = ro.r['url']
    #     gene_set_indices = limma.ids2indices(readRDS(url(gene_sets_file)), results['entrezgene_id'])
    #     return limma.camera(self.dge, gene_set_indices, self.design_df, self.contrasts, inter_gene_cor=inter_gene_cor)
    

def remove_batch(x, **kwargs):
    """Helper function to run limma's removeBatchEffect on the given data"""
    return limma.removeBatchEffect(x, **kwargs)


def cpm(dge, log=True):
    """Run edger's cpm function on given dge object"""
    return edger.cpm(dge, log=log)


class DESeq2(object):
    def __init__(self, counts_df, design_s, order, gene_name_s=None):
        """
        counts_df: counts matrix, genes x samples
        design_s:  class label for each sample
        order: test order; the first element is the reference
        """
        assert design_s.isin(order).all()
        self.method = 'DESeq2'
        self.counts_df = counts_df
        self.counts_df = counts_df[design_s.index]
        self.design_s = design_s.astype('category')
        self.design_s = self.design_s.cat.reorder_categories(order)
        self.design_s.name = 'var'
        self.design_df = pd.concat([
            pd.Series(np.ones(design_s.shape[0]), name='1', index=design_s.index),
            design_s.to_frame()
        ], axis=1)
        self.gene_name_s = gene_name_s
        self.dds = deseq2.DESeqDataSetFromMatrix(countData=self.counts_df,
            colData=ro.DataFrame(self.design_s.to_frame()),
            design=ro.Formula('~ var'))
        self.dds = deseq2.estimateSizeFactors_DESeqDataSet(self.dds)
        self.size_factors = deseq2.sizeFactors_DESeqDataSet(self.dds)

    def run(self):
        self.dds = deseq2.DESeq(self.dds)
        with localconverter(ro.default_converter + pandas2ri.converter):
            self.results_df = to_dataframe(deseq2.results(self.dds, parallel=True))

        self.results_df['baseMean'] = np.log2(self.results_df['baseMean'])
        self.results_df.rename(columns={
            'log2FoldChange':'log2FC', 'lfcSE':'log2FC_se',
            'stat':'tstat', 'baseMean':'log2_mean_expr',
            'pvalue':'pval', 'padj':'pval_adj',
        }, inplace=True)
        self.results_df.index.name = 'gene_id'
        if self.gene_name_s is not None:
            self.results_df = self.gene_name_s.rename('gene_name').to_frame().join(self.results_df, how='right')
        self.results_df = self.results_df.sort_values('pval_adj')

        self.comparison = deseq2.resultsNames(self.dds)
        print(f"DE genes at 0.05 FDR: {sum(self.results_df['pval_adj']<=0.05)}")

    def qq_plot(self):
        qtl.plot.qqplot(self.results_df['pval'], fontsize=12)

    def ma_plot(self):
        ma_plot(self.results_df['log2_mean_expr'], self.results_df['log2FC'], gridsize=100)

    def volcano_plot(self):
        volcano_plot(self.results_df['log2FC'],
                     self.results_df['pval_adj'],
                     ylabel='-log$\mathregular{_{10}}$(adj. p-value)')

    def get_counts(self, normalized=True):
        return pd.DataFrame(deseq2.counts_DESeqDataSet(self.dds, normalized=normalized),
                            index=self.counts_df.index, columns=self.counts_df.columns)


class GSEA(object):

    def __init__(self, de_results, gene_set, gsea_module='gseapy', **kwargs):
        """Object for storing GSEA prerank results and generating various plots.

        :param de_results: pandas.DataFrame with Hugo gene symbols and enrichment scores (to use for sorting)
        :param gene_set: name of gene_set collection to pass to GSEA prerank
        :param gsea_module: which GSEA library to use for computation; one of [gseapy, blitzgsea]
        
        GSEApy provides a simple interface, but is very slow for small p-values (anything with permutations > ~1000)
        blitzGSEA is very fast and accurate for small p-values, but significance estimation based on estimated gamma-distribution, not true permutation test. Also, p-values are skewed with many equal to 1 (due to weird boundary cutoffs in calculation of two-sided p-value).
        fgsea is very fast and equivalent to original GSEA code, but not implemented yet
        camera takes a different, non-parametric approach, but does allow for modeling inter-gene correlation (and is integrated well with limma voom pipeline).
        """
        self.module = gsea_module
        self.gene_set = gene_set
        self.input_ranks = de_results
        if gsea_module == 'gseapy':
            self.prerank_results = gseapy.prerank(rnk=de_results, gene_sets=gene_set, **kwargs)
            self.library = None
        elif gsea_module == 'blitzgsea':
            library = blitz.enrichr.get_library(gene_set)
            self.prerank_results = blitz.gsea(de_results, library, **kwargs)  # anchors=50, signature_cache=True, center=False
            self.library = library
        elif gsea_module == 'fgsea':
            library = blitz.enrichr.get_library(gene_set)
            self.library = library
            tmp_str = str(library).replace(': ', '=').replace('[', 'c(').replace(']', ')').replace('{', 'list(').replace('}', ')')
            library_r = ro.r(tmp_str)
            ro.r.assign('library_r', library_r)  # dict of lists to (named) list of characters
            r_de_res = ro.FloatVector(de_results['enrichment_score'].values)
            r_de_res.names = ro.StrVector(de_results['genes_hugo'].values)
            ro.r.assign('r_de_res', r_de_res)  # ???
            
            # store as an R dataframe
            self.prerank_results = fgsea.fgsea(pathways=library_r, stats=r_de_res, **kwargs)

        elif gsea_module == 'camera':
            raise(NotImplementedError('camera is not yet implemented'))

    def gseaplot(self, geneset_name, **kwargs):
        """Generate gseaplot for this geneset

        :param geneset_name: Name of geneset or index in prerank_results
        :return: Figure showing enrichment for this gene set
        """
        if self.module == 'gseapy':
            if isinstance(geneset_name, int):
                terms = self.prerank_results.res2d.Term
                geneset_name = terms[geneset_name]
            return gseapy.gseaplot(rank_metric=self.prerank_results.ranking, term=geneset_name, **self.prerank_results.results[geneset_name])
        elif self.module == 'blitzgsea':
            return blitz.plot.running_sum(self.input_ranks, geneset_name, self.library, result=self.prerank_results, **kwargs)  # compact=False, center=False

    def dotplot(self, cutoff_qval=0.2):
        """Generate dotplot for gene sets that pass the qvalue cutoff

        :param cutoff_qval: Qval cutoff (default 0.2)
        :return: Figure showing enrichment scores and q values for significant gene sets
        """
        if self.module == 'gseapy':
            return gseapy.plot.dotplot(self.get_results(), title=self.gene_set,cmap='viridis_r', size=20, figsize=(3,5), column='FDR q-val', cutoff=cutoff_qval)
        elif self.module == 'blitzgsea':
            res = self.prerank_results
            n_sig = res[res['fdr'] < cutoff_qval].shape[0]
            return blitz.plot.top_table(self.input_ranks, self.library, res, n=n_sig)

    def get_results(self):
        """Return prerank results in a pandas.DataFrame"""
        if self.module == 'gseapy':
            return self.prerank_results.res2d
        elif self.module == 'blitzgsea':
            return self.prerank_results.rename(columns={'nes': 'NES', 'es': 'ES', 'pval': 'NOM p-val', 'fdr': 'FDR q-val', 'leading_edge': 'Lead_genes'}).reset_index()
        elif self.module == 'fgsea':
            # automatic conversion breaks (use pandas to turn into dataframe)
            tmp_df_r = self.prerank_results
            non_float_df = pd.DataFrame([tmp_df_r[0], tmp_df_r[6], tmp_df_r[7]], index=[tmp_df_r.names[0], tmp_df_r.names[6], tmp_df_r.names[7]])
            float_df = pd.DataFrame(tmp_df_r[1:6], index=tmp_df_r.names[1:6], dtype=float)
            df_py = pd.concat([non_float_df.T, float_df.T], axis=1)
            return df_py.rename(columns={'pval': 'NOM p-val', 'padj': 'FDR q-val', 'leadingEdge': 'Lead_genes', 'pathway': 'Term', 'size': 'geneset_size'}).sort_values(['FDR q-val', 'NES'])
        elif self.module == 'camera':
            pass


def volcano_plot(log2_fc, pval, title='', xmax=None, ymax=None, xcut=None,
                 bc=None, ax=None, fontsize=12,
                 xlabel='log$\mathregular{_{2}}$(fold change)',
                 ylabel='-log$\mathregular{_{10}}$(p-value)'):
    """Volcano plot"""

    if ax is None:
        ax = qtl.plot.setup_figure(2,2)

    m = pval.notnull()
    ax.scatter(log2_fc[m], -np.log10(pval[m]), s=20, c='k', alpha=0.2, rasterized=True, edgecolor='none')
    qtl.plot.format_plot(ax, fontsize=fontsize-2)

    if xmax is None:
        xmax = np.ceil(np.max(np.abs(log2_fc)))
    ax.set_xlim([-xmax,xmax])
    if ymax is not None:
        ax.set_ylim([0, ymax])
    ylim = ax.get_ylim()

    if xcut is not None:  # vertical dashed lines
        ax.plot([-xcut, -xcut], ylim, '--', color=hsv_to_rgb([0,1,0.8]))
        ax.plot([xcut, xcut], ylim, '--', color=hsv_to_rgb([0,1,0.8]))

    if bc is not None:  #  Bonferroni threshold
        bc = -np.log10(0.05/len(log2_fc))
        ax.plot([-xmax, xmax], [bc,bc], '--', color=hsv_to_rgb([.6,0.6,0.8]))

    ax.spines['left'].set_position(('outward', 6))
    ax.spines['bottom'].set_position(('outward', 6))
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)


def ma_plot(log2_mu, log2_fc, highlight_ids=None, highlight_labels=None, gridsize=100):
    """MA plot"""
    cmap = copy.copy(plt.cm.get_cmap("RdYlBu_r"))
    cmap.set_bad('w', 1.)

    idx = np.isfinite(log2_fc) & np.isfinite(log2_mu)

    ax = qtl.plot.setup_figure(3, 2)
    h = ax.hexbin(log2_mu[idx], log2_fc[idx], bins='log', linewidths=1, gridsize=gridsize,
                  cmap=cmap, mincnt=1, zorder=10, edgecolors='none')#, vmin=0)
    qtl.plot.format_plot(ax, fontsize=10)
    xlim = ax.get_xlim()
    ylim = np.array(ax.get_ylim())
    b = np.max(np.abs(ylim))
    ax.plot([xlim[0]-10, xlim[1]+10], [0,0], '--', color=[0.6]*3, lw=1)

    if highlight_ids is not None:
        x = log2_mu[highlight_ids]
        y = log2_fc[highlight_ids]
        idx = np.isfinite(x) & np.isfinite(y)
        plt.scatter(x[idx], y[idx], s=24, zorder=200, color=hsv_to_rgb([0,0.9,0.9]))
        if highlight_labels is not None:
            for i,xi,yi in zip(np.array(highlight_labels)[idx], x[idx], y[idx]):
                plt.text(xi, yi, i, zorder=300, ha='right')

    ax.set_xlim(xlim)
    ax.set_ylim([-b, b])
    ax.spines['left'].set_position(('outward', 6))
    ax.spines['bottom'].set_position(('outward', 6))
    ax.set_xlabel('log$\mathregular{_2}$(mean read count)', fontsize=12)
    ax.set_ylabel('log$\mathregular{_2}$(fold change)', fontsize=12)


def compare_effects(de1, de2):
    """Compare estimated effect sizes (logFC)"""
    ix = de1.results_df.index
    ax, cax = qtl.plot.hexdensity(de1.results_df.loc[ix, 'log2FC'],
                                  de2.results_df.loc[ix, 'log2FC'],
                                  scale='linear', gridsize=100)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    lims = [np.minimum(xlim[0], ylim[0]), np.maximum(xlim[1], ylim[1])]
    ax.plot(lims, lims, 'k--', zorder=0, alpha=0.33)

    ax.set_xlabel(de1.method+' [-log$\mathregular{_2}$(FC)]', fontsize=12)
    ax.set_ylabel(de2.method+' [-log$\mathregular{_2}$(FC)]', fontsize=12)
