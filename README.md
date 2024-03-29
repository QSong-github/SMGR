# SMGR: a joint statistical method for integrative analysis of single-cell multi-omics data

[![DOI](https://zenodo.org/badge/405528473.svg)](https://zenodo.org/badge/latestdoi/405528473)

A tailored tool for identifying coherent signals across single-cell multi-omics data. 

Unravelling the regulatory programs from single-cell multi-omics data has been one of the major challenges, especially in the current emerging single-cell field. Currently there is a gap between fast-growing single-cell multi-omics data and effective methods for the integrative analysis of these inherent sparse and heterogeneous data. Here we proposed SMGR to detect coherent functional regulatory signals and target genes from the joint scRNA-seq and scATAC-seq data obtained from different samples. SMGR enables the identification of co-regulatory programs and the elucidation of regulating mechanisms.

## Installation
```
devtools::install_github("QSong-github/SMGR")
```

## How to use

SMGR tutorial provides examples and explanations of its functions and how to use them. This documentation introduces the main features of SMGR.
- [SMGR tutorial](https://github.com/QSong-github/SMGR/blob/main/vignette/SMGR_vignettes.pdf)

## load SMGR package
```
library(SMGR)
```
## input data: a list of scRNA-seq data and scATAC-seq data

SMGR works with a list of scRNA-seq and scATAC-seq dataset as inputs. Bascially, the format of data input is as follows. Example data files can be found in the ```data``` folder.

For scRNA-seq data:

| CellID | rna-cell1 | rna-cell2 | rna-cell3 | rna-cell4 | ... |
|----|--------|--------|--------|---------|-----|
| Gene1 | 0 | 1 | 0 | 0 | ... |
| Gene2 | 0 | 0 | 1 | 0 | ... |
| Gene3 | 0 | 0| 1 | 0  | ... |
|...    |...|...|...|...|...|

For scATAC-seq data:

| CellID | atac-cell1 | atac-cell2 | atac-cell3 | atac-cell4 | ... |
|----|--------|--------|--------|---------|-----|
| Peak1 | 1 | 0 | 0 | 0 | ... |
| Peak2 | 1 | 0 | 1 | 0 | ... |
| Peak3 | 0 | 0| 0 | 0  | ... |
|...    |...|...|...|...|...|

Load the example data
```
rna.cts <- readRDS('./data/simuation_scRNA-seq.RDS')
atac.cts <- readRDS(file='./data/simuation_scATAC-seq.RDS')
input_data <- list(rna.cts, atac.cts)
```
## SMGR process

Input data is a list of scRNA-seq and scATAC-seq data

```
result1 <- smgr_main(sm.data = input_data, K=nrow(input_data[[1]]))
```
result1 contains the latent representation of joint scRNA-seq and scATAC-seq data

For details, please refer to [SMGR tutorial](https://github.com/QSong-github/SMGR/blob/main/vignette/SMGR_vignettes.pdf).

## Cite

Please cite our paper if you use this code in your own work:

```
Song Q, Zhu X, Jin L, Chen M, Zhang W, Su J. "SMGR: a joint statistical method for integrative analysis of single-cell multi-omics data." NAR genomics and bioinformatics. 2022 Sep 1;4(3):lqac056.
```
