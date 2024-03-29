
#' @title SMGR main function
#' @description extract co-regulation genes across single-cell multi-omics data
#' @param sm.data consists of single-cell RNA-seq data and single-cell ATAC-seq data in the format of matrix. Columns represent cells and rows represent genes (scATAC-seq input has been converted to gene-based data matrix)
#' @param K The number of features of the input data
#' @param R the dimension of latent representation 
#' @param n_burnin burnin for mcmc sampling
#' @param n_draw drawing for mcmc sampling
#' @param n_iter number of iterations
#' @param option nb or zinb
#' @param type simulate or experimental data
#' @return latent representation; co-regulation programs
#' @export smgr_main function

smgr_main <- function(sm.data, K, R = 3, n_burnin=200, n_draw=200, n_iter=20, option='nb', eps=1e-4,
                      type='simulate', output = 'latent')
{
    pkg = c('glmnet','MASS','purrr','mpath','zic','pscl','parallel')
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)>0){
      BiocManager::install(new.pkg, dependencies = TRUE)
    } else {
      sapply(pkg, require, character.only = TRUE)
    }
  
    cat("## ============================================================================\n")
    cat("## |                             running SMGR                                  |\n")
    cat("## ----------------------------------------------------------------------------\n")

    len <- length(sm.data)
    
    sm.datas <- transform(sm.data, len, R, type) 
    
    if (len==2){
        cat("## |            scRNA-seq data and scATAC-seq data identified ....              |\n")
    } else if (len<2){
          cat("## |               single-cell multi-omics data is not identified ....       |\n")}

    cat("## |                           Set initial values ....                         |\n")
    mean_X = matrix(0, nrow = K, ncol = R); dif = 0.1
    init_Y <- matrix(rnorm(K*R,0,1), nrow = K, ncol = R)
    update_Y <- init_Y

    iter = 0

    while (iter < n_iter && dif> eps)
    {
        iter = iter + 1
        
  if (iter<10){term <- paste0('iter= ',iter)}else {term <- paste0('iter=',iter)}
           cat(paste0('## |                               ',term,'.....                                |\n'))

        if (option=='zinb'){

            estimates <- mclapply(sm.datas, function(dat){
                beta0 <- dat$b0; beta1 <- dat$b1

                estimate_i <- mclapply(1:ncol(dat$x_i), function(f)
                {
                    datas <- data.frame(update_Y,y=dat$x_i[,f]);

                    tmp <- mclapply(1:(ncol(datas)-1),function(l){
                        m1 <- zeroinfl(y ~ ., dist = 'negbin',
                                       data = datas[,c(l,ncol(datas),drop=FALSE)])
                        return(m1)}, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)

                    incept1 <- mean(sapply(tmp,function(m1){coef(m1)[1]}))
                    incept2 <- as.numeric(sapply(tmp,function(m1){coef(m1)[2]}))
                    zinb.coef <- as.matrix(c(incept1,incept2))
                    zinb.th <- mean(sapply(tmp,function(m1){summary(m1)$theta}))
                    zinb.mu <- as.numeric(sapply(tmp,function(m1){m1$fitted.values[which.min(m1$residuals)]}))
                    temp <- list(zinb_coef = zinb.coef,
                                 zinb_th = zinb.th,
                                 zinb_mu = zinb.mu)
                    return(temp)
                }, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)
                
                dat$coef <- do.call(cbind, map(estimate_i, 1))
                dat$th <- do.call(cbind, map(estimate_i, 2))
                dat$mu <- do.call(cbind, map(estimate_i, 3))
                dat$b0 <- dat$coef[1,]; dat$b1 <- t(dat$coef[-1,])
                dat$dif = abs(cbind(dx$b0-beta0,dx$b1-beta1))
                return (dat)
            }, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)
            }
        estimates <- mclapply(sm.datas, function(dat){
            beta0 <- dat$b0; beta1 <- dat$b1
            estimate_i <- mclapply(1:ncol(dat$x_i), function(f)
                {
                    datas <- data.frame(update_Y,y=dat$x_i[,f]);
                    nb_reg <- glmregNB(y ~ ., data=datas, family='negbin', penalty="enet",
                                       lambda =1, alpha=1, thresh=1e-4,
                                       standardize=FALSE, parallel=FALSE);

                    temp <- list(nb_coef = as.matrix(coef(nb_reg)),
                                 nb_th = nb_reg$theta,
                                 nb_mu = nb_reg$fitted.values)
                    return(temp)}, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)

            dat$coef <- do.call(cbind, map(estimate_i, 1))
            dat$th <- do.call(cbind, map(estimate_i, 2))
            dat$mu <- do.call(cbind, map(estimate_i, 3))
            dat$b0 <- dat$coef[1,]; dat$b1 <- t(dat$coef[-1,])
            return (dat)
        }, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)
        
        tmp_update <- tryCatch(
        {
            update_Ys(mean_x = mean_X,
                      update_y = update_Y,
                      estimate_y = estimates,
                      k = K,
                      r = R,
                      n = len,
                      nburnin = n_burnin,
                      ndraw = n_draw)
        },
        error=function(cond) {
            message("\n ##  Errors occured during optimization... \n")
            message("\n ## Here's the original error message... \n")
            message(cond)
            return(NA)
        },
        warning=function(cond) {
                return(NULL)
                }
        )

        update_Y <- tmp_update
        }


    kmeans.fit = kmeans(update_Y, R, nstart=100)
    Yclusters = kmeans.fit$cluster

    if (output =='all'){
        totalBic <- sum(unlist(
            mclapply(estimates, function(dx)
        {
            BICs(xi = dx$x_i, beta = dx$b1, th = dx$th[1,], mu = dx$mu, update = update_Y)
        }, mc.cores = detectCores(), mc.preschedule = FALSE, mc.set.seed = FALSE)));
    latents <- list(estimate_f = estimates, latent_f = update_Y, bic_f = totalBic)
    } else {
        latents <- list(latent_f = update_Y, clusters = Yclusters) }
        
    return (latents)

    cat("## ======================= finished ===================================\n")
}

