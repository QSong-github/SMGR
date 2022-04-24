# packages = c('glmnet','MASS','purrr','mpath','zic','pscl','parallel')
# # BiocManager::install(packages)
# 
# # devtools::install_github("QSong-github/SMGR")
# library('SMGR')
# 
# # Note: please invoke these packages to run the following codes:
# 
# # Run SMGR with example data
# 
# # Example data is deposited in the data folder.
# # @param rna.cts simulated scRNA-seq data
# rna.cts <- readRDS('./data/simuation_scRNA-seq.RDS')
# 
# # @param atac.cts simulated scATAC-seq data
# atac.cts <- readRDS(file='./data/simuation_scATAC-seq.RDS')
# 
# input_data <- list(rna.cts, atac.cts)
# 
# #**SMGR process**
# 
# # Input data is a list of scRNA-seq and scATAC-seq data
# input_data <- list(as.matrix(rna.cts),as.matrix(atac.cts))
# 
# result1 <- smgr_main(sm.data = input_data, K=nrow(input_data[[1]]))
# 
# # member1: is the ground truth
# member1 <- c(rep(1,600),rep(2,600),rep(3,600))
# 
# # programs: co-expressed gene & peaks from latent representation
# programs <- result1$clusters
