# SMGR

A tailored tool for identifying coherent signals and co-regulators across single-cell multi-omics data

# 1. Installation
```
devtools::install_github("QSong-github/SMGR")
```

# 2. How to use

## 2.1 load SMGR package
```
library(SMGR)

```
SMGR works with single-cell RNA-seq dataset and single-cell ATAC-seq datasets as inputs. Example data is shown in ```Data``` folder.

## 2.2 input data list with scRNA-seq and scATAC-seq datasets

SMGR works with a list of scRNA-seq and scATAC-seq dataset as inputs. Bascially, the format of data input is as follows. Example data files can be found in the ```Data``` folder.

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

#' load the example data
```
data("data1")
```
data1 is the data list of scRNA-seq and scATAC-seq data

## 2.3 run main function of the example data
```
result1 <- smgr_main(sm.data = data1, K=3, N=nrow(data1[[1]]))
```
result1 contains the latent representation of joint scRNA-seq and scATAC-seq data

## 2.4 evaluate of clustering results using ground truth (this is optional)

calculate the Adjusted Rand Index
```
library(clues)
adjustedRand(result1$clusters,example1.member)
```
## 3. Examples and reproducible results 

Examples can be found using the example.R script

Identify the optimal co-regulation programs with least BIC value
```
files <- list.files(path=path1, pattern='results.RDS')

opt <- vector('list')
for ( i in 1:length(files)){
    opt[[i]] <- readRDS(files[i])
    names(opt[[count]]) <- paste0('opt_',strsplit(strsplit(i,'_')[[1]][4],'results')[[1]][1])}
```

Identify the result with least BIC value
```
bics <- sapply(1:length(opt),function(i){ f.opt <- opt[[i]][[1]]$BIC })
optimal <- grep(min(bics),bics)
```
