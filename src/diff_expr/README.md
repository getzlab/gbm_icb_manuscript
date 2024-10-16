## Differential expression testing

This repository contains a Python module with wrapper classes for the following methods (R packages):
* [limma-voom](https://bioconductor.org/packages/release/bioc/html/limma.html); [Law et al., Genome Biol, 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html); [Love et al., Genome Biol, 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)


### Requirements

This module requires the following Python modules: `pandas`, `numpy`, `matplotlib`, `rpy2`, `qtl`, `gseapy`

An R installation with the following packages is also needed: `limma`, `edgeR`, `DESeq2`

### Installation

After cloning this package (`git clone git@github.com:ConorMesser/diff_expr.git`), create a conda environment to run analyses in using the following commands:

```
conda create --name <env_name> -c conda-forge -c bioconda python=3.8 r-base=4.1.1 bioconductor-limma bioconductor-edger bioconductor-deseq2 r-smartsva r-isva gseapy pandas numpy matplotlib
conda activate <env_name>
pip install qtl rpy2
```

These packages can also be installed to an existing conda environment, however dependency errors may arise.

Alternatively, one can create a new conda environment using the provided environment.yml file:

```conda env create --name <env_name> -f ./environment.yml```