# SMGR

A tailored tool for identifying coherent signals and co-regulators across single-cell multi-omics data

## Installation
```
devtools::install_github("QSong-github/SMGR")
```

## How to use

SMGR works with single-cell RNA-seq dataset and single-cell ATAC-seq datasets as inputs. Example data is shown in ```Data``` folder.

SMGR tutorial provides examples and explanations of its functions and how to use them. This documentation introduces the main features of SMGR.
* [SMGR tutorial](https://github.com/QSong-github/SMGR/blob/main/vignette/vignette_make.html)
*
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

## For details, please refer to * [SMGR tutorial](https://github.com/QSong-github/SMGR/blob/main/vignette/vignette_make.html)
