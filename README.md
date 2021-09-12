# SMGR

A tailored tool for identifying coherent signals and co-regulators across single-cell multi-omics data

# 1. Installation
```
devtools::install_github("QSong-github/SMGR")
```

# 3. How to use

## 3.1 load SMGR package
```
library(SMGR)

```
SMGR works with single-cell RNA-seq dataset and single-cell ATAC-seq datasets as inputs. Example data is shown in ```Data``` folder.

## 3.2 input data list with scRNA-seq and scATAC-seq datasets
```
#' load the example data
data("data1")

```

## 3.3 run main function of the example data
```
result1 <- smgr_main(sm.data = data1, K=3, N=nrow(data1[[1]]))

```
result1 contains the latent representation of joint scRNA-seq and scATAC-seq data
