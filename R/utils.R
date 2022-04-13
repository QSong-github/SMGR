
#' @title calculate the Dunn index and Sihouette index
#' @param mat input data
#' @param cluster identified programs
#' @export evaluation metrics
Evaluate <- function(mat,cluster){
    if (sum(is.na(rownames(mat)))>0){
        ps <- which(is.na(rownames(mat)))
        mat <- mat[-ps,]
        cluster <- cluster[-ps]}
    Dist1 <- dist(mat, method = "euclidean")
    Dist2 <- as.dist(1 - cor(t(mat)))
    out <- t(data.frame(Dunn = dunn(Dist1,cluster),
                        Sihouette = mean(silhouette(cluster, dmatrix=as.matrix(Dist1))[,3])))
                                                    
    return(out)}


#' @title Gibbs sampling for identifying latent representation
#' @param mean_x A non-negative integer. or vector
#' @return updated latent representation
#' @export  update_Ys function
update_Ys <- function(mean_x,
                      update_y,
                      estimate_y,
                      k,
                      r,
                      n,
                      sdev=0.05,
                      nburnin,
                      ndraw) {

    a = as.list(1:2)
    b = as.list(1:2)
    x <- lapply(1:2,function(i){matrix(0, nrow = k, ncol = r)})
    mu <- lapply(1:2,function(i){matrix(0, nrow = k, ncol = r)})
    th <- lapply(1:2,function(i){matrix(0, nrow = k, ncol = r)})
    p = rep(1:2)

    for ( j in 1:n){
        a[[j]] = estimate_y[[j]]$b0
        b[[j]] = estimate_y[[j]]$b1
        x[[j]] = estimate_y[[j]]$x_i
        mu[[j]] = estimate_y[[j]]$mu
        th[[j]] = estimate_y[[j]]$th
        p[j] = estimate_y[[j]]$p
    }

  temp <- do.call(rbind,
                  mclapply(1:k,function(t)
                  {
                      .C("gibbs", meanz = as.double(mean_x[t,]),
                         lastz = as.double(update_y[t,]),
                         as.integer(nburnin),
                         as.integer(ndraw),
                         as.double(a[[1]]),
                         as.double(b[[1]]),
                         as.integer(x[[1]][t,]),
                         as.double(mu[[1]][t,]),
                         as.double(th[[1]][1,]),
                         as.integer(p[1]),
                         as.double(a[[2]]),
                         as.double(b[[2]]),
                         as.integer(x[[2]][t,]),
                         as.double(mu[[2]][t,]),
                         as.double(th[[2]][1,]),
                         as.integer(p[2]),
                         as.integer(r),
                         as.double(sdev),
                         as.integer(n),PACKAGE = 'SMGR')[[1]]
                  }, mc.cores = detectCores(), mc.preschedule = FALSE, mc.set.seed = FALSE))

  return (temp)}

#' @title calculate the BIc value
#' @description calculate bic
#' @details calculate the BIc value to seletct the best latent variables
#' @param xi
#' @return the BIC value
#' @export  BICs function

BICs <- function(xi, beta, th, mu, update){
  BIC = 0
  all.ll <- sum(unlist(lapply(1:nrow(xi),function(i){
    loglik(th,mu[i,],y=xi[i,])})))
  BIC = BIC - 2*all.ll + (sum(beta!=0))*log(nrow(xi)*ncol(xi))
  return (BIC)}

#' @title calculate the loglikehood
#' @description loglikehood
#' @details calculate the loglikehood based on the negative binomial distribution
#' @param th theta in the negative binomial distribution
#' @param mu mu in the negative binomial distribution
#' @return A numeric vector of log density.
#' @export  posterior likelihood

loglik <- function(th, mu, y){
  sum(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
        y * log(mu + (y == 0)) - (th + y) * log(th + mu))
}
