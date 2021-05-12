mOrdinaryBootCV <- function(n, p, MaxLag, alpha, b, parallel = FALSE) {

  x <- Rfast::matrnorm(n, p)
  if ( missing(MaxLag) || MaxLag < 0 )  stop( "'MaxLag' must be greater than 0." )

  mADCFfun <- function(x, lags) {
    mat <- matrix(NA, nrow = p, ncol = p)
    X <- x[(1 + lags):n, , drop = FALSE]
    Y <- x[1:(n - lags), , drop = FALSE]
    for ( i in 1:p )  mat[i, ] <- dcov::mdcor(X[, i], Y)
    return( sqrt(mat) )
  }

  rstar <- function(k){
    cv <- vector()
    Rstark <- function(k) {
      ind <- sample(1:n, replace = TRUE)
      xStar <- x[ind, ]
      mADCFfun(xStar, lags = k)
    }
    result <- replicate(b, Rstark(k) )

    ela <- expand.grid(1:p, 1:p)
    for ( i in 1:dim(ela)[1] ) {
      quant <- quantile( result[ela[i, 1], ela[i, 2], ], 1 - alpha )
      pv <- mean( result[ela[i, 1], ela[i, 2], ] >= quant )
      pvadj <- p.adjust( pv, method = "fdr" )
      cv[i] <- quantile( result[ela[i, 1], ela[i, 2], ], 1 - pvadj )
    }
    cv
  }

  if ( parallel ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    clusterSetRNGStream(cl = cl, iseed = 9182)
    i <- 1:MaxLag
    fe_call <- as.call( c( list( as.name("foreach"), i = i, .combine="rbind",
                          .export = "mADCF", .packages = "dcov") ) )
    fe <- eval(fe_call)
    Rstar <- fe %dopar% rstar(i)
    stopCluster(cl)
    mcv <- matrix( sapply(1:p^2, function(i) max( Rstar[, i] ) ), nrow = p, ncol = p, byrow = TRUE)
  } else {
    res <- lapply(1:MaxLag, FUN = function(i) rstar(i) )
    Rstar <- t( sapply(1:2, FUN = function(i) unlist( res[[ i ]]) ) )
    mcv <- matrix( sapply(1:p^2, function(i) max( Rstar[[ i ]] ) ), nrow = p, ncol = p, byrow = TRUE )
  }

  return(mcv)
}














# mOrdinaryBootCV <- function(n,q,MaxLag,b,parallel=FALSE){
#  x <- replicate(q,rnorm(n))
#  if (missing(MaxLag) || MaxLag < 0)
#       stop("'MaxLag' must be greater than 0")
#  mADCFfun <- function(x,lags){
#    mat <- matrix(NA,nrow=q,ncol=q)
#    X <- rbind(x[(1+lags):n,])
#    Y <- rbind(x[1:(n-lags),])
#   for (i in 1:q){
#    for (j in 1:q){
#     mat[i,j] <- dcor(X[,i],Y[,j])
#    }
#   }
#   return(mat)
#  }
#   rstar <- function(k){
#    cv <- vector()
#    Rstark <- function(k){
#     ind <- sample(1:n,replace=T)
#     xStar <- x[ind,]
#     Rm <- mADCFfun(xStar,lags=k)  #### wrong here, it should be xStat, not x
#     return(Rm)
#    }
#    result <- replicate(b,Rstark(k))
#    s <- 1
#    for(i in 1:q){
#     for (j in 1:q){
#      quant <- quantile(result[i,j,],0.95)
#      pv <- mean(result[i,j,]>=quant)
#      pvadj <- p.adjust(pv,method="fdr")
#      cv[s] <- quantile(result[i,j,],1-pvadj)
#      s <- s+1
#     }
#    }
#    return(cv)
#   }
#  if(parallel==TRUE){
#   closeAllConnections()
#   cl <- makeCluster(2)
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- 1:MaxLag
#   fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="rbind",.export=c("mADCF","dcor")) ))
#   fe <- eval(fe_call)
#   Rstar <- fe %dopar% rstar(i)
#   stopCluster(cl)
#   mcv <- matrix(sapply(1:q^2, function(i) max(Rstar[,i])),nrow=q,ncol=q,byrow=T)
#  }
#  else {
#  res <- lapply(1:MaxLag,FUN=function(i) rstar(i))
#  Rstar <- t(sapply(1:2,FUN=function(i) unlist(res[[i]])))
#  mcv <- matrix(sapply(1:q^2, function(i) max(Rstar[[i]])),nrow=q,ncol=q,byrow=T)
#  }
#  return(mcv)
# }
