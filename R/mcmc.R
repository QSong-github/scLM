
#' @title convert data
#' @param datatype
#' @export  data conversion

Convert <- function(datatype,ndt,K){
  if (ndt>0){
    dat <- lapply(datatype,function(x){
      meanx = apply(x,2,mean);
      B0 = log(meanx); B1=matrix(0,nrow=ncol(x),ncol=K)
      p=ncol(x)
      list(xi=x,meanx=meanx,B0=B0,B1=B1,p=p)})
    return (dat)} else{
      stop("Error: at least one dataset")}}

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

#' @title MCMC sampling
#' @description  update the latent variables
#' @param x A non-negative integer. or vector
#' @return updated latent variables values
#' @export  updateZ function

updateZ <- function(meanX,finalZ,Res,N,K,ndt,sdev,n.burnin,n.draw){
  a = as.list(1:9)
  b = as.list(1:9)
  x <- lapply(1:9,function(i){matrix(0,nrow=N,ncol=K)})
  mu <- lapply(1:9,function(i){matrix(0,nrow=N,ncol=K)})
  th <- lapply(1:9,function(i){matrix(0,nrow=1,ncol=K)})
  p = rep(1:9)

  for ( o in 1:ndt){
    a[[o]] = Res[[o]]$beta0
    b[[o]] = Res[[o]]$beta1
    x[[o]] = Res[[o]]$xi
    mu[[o]] = Res[[o]]$Mus
    th[[o]] = Res[[o]]$Ths
    p[o] = Res[[o]]$p
  }

  temp <- do.call(rbind,
                  mclapply(1:N,function(t)
                  {
                    .C("nbMcmc",meanz=as.double(meanX[t,]),lastz=as.double(finalZ[t,]),as.integer(n.burnin),as.integer(n.draw),
                       as.double(a[[1]]),as.double(b[[1]]),as.integer(x[[1]][t,]),as.double(mu[[1]][t,]),as.double(th[[1]][1,]),as.integer(p[1]),
                       as.double(a[[2]]),as.double(b[[2]]),as.integer(x[[2]][t,]),as.double(mu[[2]][t,]),as.double(th[[2]][1,]),as.integer(p[2]),
                       as.double(a[[3]]),as.double(b[[3]]),as.integer(x[[3]][t,]),as.double(mu[[3]][t,]),as.double(th[[3]][1,]),as.integer(p[3]),
                       as.double(a[[4]]),as.double(b[[4]]),as.integer(x[[4]][t,]),as.double(mu[[4]][t,]),as.double(th[[4]][1,]),as.integer(p[4]),
                       as.double(a[[5]]),as.double(b[[5]]),as.integer(x[[5]][t,]),as.double(mu[[5]][t,]),as.double(th[[5]][1,]),as.integer(p[5]),
                       as.double(a[[6]]),as.double(b[[6]]),as.integer(x[[6]][t,]),as.double(mu[[6]][t,]),as.double(th[[6]][1,]),as.integer(p[6]),
                       as.double(a[[7]]),as.double(b[[7]]),as.integer(x[[7]][t,]),as.double(mu[[7]][t,]),as.double(th[[7]][1,]),as.integer(p[7]),
                       as.double(a[[8]]),as.double(b[[8]]),as.integer(x[[8]][t,]),as.double(mu[[8]][t,]),as.double(th[[8]][1,]),as.integer(p[8]),
                       as.double(a[[9]]),as.double(b[[9]]),as.integer(x[[9]][t,]),as.double(mu[[9]][t,]),as.double(th[[9]][1,]),as.integer(p[9]),
                       as.integer(K),as.double(sdev),as.integer(ndt),PACKAGE = 'scLM')[[1]] }, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE))

  return (temp)}

#' @title calculate the BIc value
#' @description calculate bic
#' @details calculate the BIc value to seletct the best latent variables
#' @param xi
#' @return the BIC value
#' @export  BICs function

BICs <- function(xi,beta,th,mu,finalZ){
  BIC = 0
  all.ll <- sum(unlist(mclapply(1:nrow(xi),function(i){
    loglik(th,mu[i,],y=xi[i,])}, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)))
  BIC = BIC - 2*all.ll + (sum(beta!=0))*log(nrow(xi)*ncol(xi))
  return (BIC)}
