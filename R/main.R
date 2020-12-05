
#' required packages
#' @import from utility.R
library(mpath); library(zic); library(pscl); library(glmnet); library(MASS); library(parallel); library(purrr)

#' @title scLM main function
#' @description extract co-expressed genes across multiple single-cell datasets simultaneously
#' @param datalist consists of multiple scRNA-seq datasets in the format of matrix.Columns correpond to objects (cell barcodes for example), rows correspond to attributes or variables (gene symbols for example).
#' @param N The number of rows in the data matrix
#' @param K The number of co-expression clusters in the single cell data
#' @param n.burnin the number of burnin for the MCMC process
#' @param n.draw the number of drawing for the MCMC process
#' @param maxiter An integer > 0. Default is 20. The number of iterations.
#' Outputs:
#' @return the latent variables; coefficients; co-expression clusters; BIC values
#' @export Multi_NB function
#' @example
#' \dontrun{
#' Pass a value to node as a replacement for FindAllMarkersNode
#' system.file("extdata", "2012.csv", package = "testdat")}
#' Multi_NB(datalist=trial, N=nrow(trial[[1]]),K=3)


Multi_NB <- function(datalist,N,K,n.burnin=200,n.draw=200,maxiter=20,eps=1.0e-4,sdev=0.05,choice=NULL,lambda=1,alpha=1,thresh=1e-04)
{

    ndt <- length(datalist)
    Datas <- Convert(datalist,ndt,K)

    maxdif <- 1
    iter = 0;

    meanX = matrix(0,nrow=N,ncol=K)
    initZ <- matrix(rnorm(N*K,0,1),nrow=N,ncol=K)
    finalZ <- initZ

    cat("## ============================================================================\n")
    cat("## |                             PROCESSING                                    |\n")
    cat("## ----------------------------------------------------------------------------\n")
    cat("## |                           Set initial values ....                         |\n")

    if (ndt==1){
        cat("## |                      One sample identified ....                           |\n")
    } else if (ndt>1){
          cat("## |                       Multiple samples identified ....                    |\n")}

    if (is.null(choice)) #' by default lambda=1, which is faster
        {
    cat("## |                          Choose default lambda ....                       |\n")}

    while(iter<maxiter && maxdif> eps)
        {
        iter = iter + 1

        #' cross-validation for best labmda
        if(!is.null(choice)){

            cat("\n ## seeking for best labmda .... \n")

        Res <- mclapply(Datas,function(dx){
                   beta0 <- dx$B0; beta1 <- dx$B1
                   res <- mclapply(1:ncol(dx$xi),function(f)
                       {
    datas <- data.frame(finalZ,y=dx$xi[,f]);
    fm_nb <- glmregNB(y ~ .,data=datas,family='negbin',penalty="enet",alpha=1,thresh=1e-4,standardize=FALSE);
    minBic <- which.min(BIC(fm_nb));
    return(list(coefs=as.matrix(coef(fm_nb, minBic)),ths=fm_nb$theta[minBic],mus=fm_nb$fitted.values[,minBic],lambdas=fm_nb$lambda[minBic]))},mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)

                   dx$Coef <- do.call(cbind,map(res, 1))
                   dx$Ths <- do.call(cbind,map(res, 2))
                   dx$Mus <- do.call(cbind,map(res, 3))
                   dx$beta0 <- dx$Coef[1,]; dx$beta1 <- t(dx$Coef[-1,])
                   dx$dif <- abs(cbind(dx$beta0-beta0,dx$beta1-beta1))
                   dx$B0 <- dx$beta0; dx$B1 <- dx$beta1
                   dx$maxdif <- max(as.vector(dx$dif))
                       dx$lambda.mins <- do.call(cbind,map(res, 4))
                   return (dx) }, mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)

    } else if (is.null(choice)) #' by default lambda=1, which is faster
          {
    if (iter<10){term <- paste0('iter= ',iter)}else {term <- paste0('iter=',iter)}
    cat(paste0('## |                               ',term,'.....                                |\n'))
            
    library(future.apply); plan(multisession);    
    Res <- future_lapply(Datas,function(dx){
                        beta0 <- dx$B0; beta1 <- dx$B1
                        res = future_lapply(1:ncol(dx$xi),FUN=function(f)
                            {
    datas <- data.frame(finalZ,y=dx$xi[,f]);
    fm_n1 <- glmregNB(y ~ .,data=datas,family='negbin',penalty="enet",lambda=lambda,alpha=alpha,thresh=thresh,standardize=FALSE);
    return(list(coefs=as.matrix(coef(fm_n1)),ths=fm_n1$theta,mus=fm_n1$fitted.values))})

                        dx$Coef <- do.call(cbind,map(res, 1))
                        dx$Ths <- do.call(cbind,map(res, 2))
                        dx$Mus <- do.call(cbind,map(res, 3))
                        dx$beta0 <- dx$Coef[1,]; dx$beta1 <- t(dx$Coef[-1,])
                        dx$dif <- abs(cbind(dx$beta0-beta0,dx$beta1-beta1))
                        dx$B0 <- dx$beta0; dx$B1 <- dx$beta1
                        dx$maxdif <- max(as.vector(dx$dif))
                        dx$lambda.mins <- 1
                        return (dx) })
}

        out <- tryCatch({updateZ(meanX,finalZ,Res,N,K,ndt,sdev,n.burnin,n.draw)},
                           error=function(cond) {
                               message("\n ##  Errors occured in the MCMC process... \n")
                               message("\n ## Here's the original error message... \n")
                               message(cond)
                               #' return NA in case of warning
                               return(NA)
                           },
                           warning=function(cond) {
                               message("\n ## Warning occured in the MCMC process... \n")
                               message("\n ## Here's the original Warning message... \n")
                               message(cond)
                               #' return NULL in case of warning
                               return(NULL)})

        if (is.na(out)){ stop("Errors occur in MCMC")} else { finalZ <- out }                                          }

    kmeans.fit=kmeans(finalZ,K,nstart=100)
    clusters=kmeans.fit$cluster
    centers=kmeans.fit$centers

    totalBic <- sum(unlist(mclapply(Res,function(dx){
                                        BICs(xi=dx$xi,beta=dx$beta1,th=dx$Ths[1,],mu=dx$Mus,finalZ=finalZ)},
                                    mc.cores = detectCores(), mc.preschedule=FALSE,mc.set.seed=FALSE)))
    return (list(DataRes=Res, clusters=clusters, centers=centers, meanZ=finalZ, BIC=totalBic))
    cat("## ============================================================================\n")
}

