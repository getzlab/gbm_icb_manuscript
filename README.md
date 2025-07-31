# gbm\_icb\_manuscript

Code for generating figures in the Glioblastoma ICB manuscript of Ghannam et al.:

> "Tumor Transcriptional Subtype is a Major Axis of Immune Checkpoint Blockade Response in Glioblastoma". _Ghannam, Iorgulescu, Weiss, Merrell, Wu et al.(2024)_ 


# Regenerating the figures:

1. Clone this Git repository
2. Install python dependencies (listed in `environment.yaml`)
3. Download data from Zenodo (DOI: 10.5281/zenodo.13942487)
4. Unzip the data in this directory: `tar -xvzf gbm_icb_data.tar.gz`. At this point there should be a directory called `data` in the `ghannam_gbm_icb_manuscript` directory.
5. Run each of the Jupyter (`*.ipynb`) and Rmarkdown (`*.Rmd`) notebooks in `notebooks/`.
   Each notebook generates one or more of the figures.
   (The `src` directory contains some helper functions employed by the Jupyer notebooks.)
