RbootCV <- function(n, MaxLag, alpha = 0.05, b, parallel = FALSE) {

  if ( missing(MaxLag) || MaxLag < 0 )  stop( "'MaxLag' must be greater than 1." )

  if (n < 500) {
    x <- rnorm(n)
  } else x <- Rfast::Rnorm(n)
  A0 <- Rfast::vecdist(x)
  Atilde0 <- ATilde(A0)
  dvarx <- mean(Atilde0 * Atilde0)

  rstar <- function(k, dvarx) {
    if ( b %% 2 == 0 ) {
      b1 <- b2 <- b/2
    } else {
      b1 <- floor(b/2)
      b2 <- ceiling(b/2)
    }
    cross <- crossDist(x, lags = k)
    Atilde <- ATilde( cross$A[[ 1 ]] )
    Btilde <- ATilde( cross$B[[ 1 ]] )
    com <- Atilde * Btilde
    Wtstar <- Rfast::matrnorm(b1, n - k) ## matrix of Z variables with b/2 rows and n - k columns
    V <- Rfast::rowsums( Wtstar %*% com * Wtstar ) / dvarx  ## this is in fact a Mahalanobis distance
    Wtstar <- Rfast::matrnorm(b2, n - k)
    V2 <- Rfast::rowsums( Wtstar %*% com * Wtstar ) / dvarx
    V <- c(V, V2)
    return( sqrt(V) / (n - k) )
  }

  if ( parallel ) {
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    clusterSetRNGStream(cl = cl, iseed = 9182)
    i <- 1:MaxLag
    fe_call <- as.call( c( list( as.name("foreach"), i = i, .export = c("crossDist", "ATilde"), .packages = "Rfast" ) ) )
    fe <- eval(fe_call)
    Rstar <- fe %dopar% rstar(i, dvarx = dvarx)
    stopCluster(cl)
    Rstar <- matrix( unlist(Rstar), ncol = MaxLag )

  } else  {
    Rstar <- matrix(nrow = b, ncol = MaxLag)
    for ( j in 1:MaxLag )  Rstar[, j] <- rstar(j, dvarx = dvarx)
  }

  quant <- Rfast2::colQuantile(Rstar, probs = 1 - alpha)
  pv <- pvadj <- Rfast::colmeans( Rfast::eachrow( Rstar, quant, oper = "-") >= 0)
  for (i in 1:MaxLag)  pvadj[i] <- p.adjust(pv[i], method = "fdr")
  cvadj <- Rfast2::colQuantile(Rstar, 1 - pvadj)

  max(cvadj)
}










# RbootCV <- function(n,MaxLag,b,parallel=FALSE){
#  x <- rnorm(n)
#  if (missing(MaxLag) || MaxLag < 0)
#       stop("'MaxLag' must be greater than 1")
#  A0 <- crossDist(x,lags=0)$A[[1]]
#  Atilde0 <- ATilde(A0)
#  rstar <- function(k){
#   cross <- crossDist(x,lags=k)
#   A <- cross$A[[1]]
#   B <- cross$B[[1]]
#   Atilde <- ATilde(A)
#   Btilde <- ATilde(B)
#   Rstark = function(Atilde,Btilde,k){
#    Wtstar <- rbind(rnorm(n-k))
#    dcov <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))
#    dvarx <- sqrt(mean((Atilde0*Atilde0))*mean((Atilde0*Atilde0)))
#    return(dcov/sqrt(dvarx))
#   }
#  return(replicate(b,Rstark(Atilde,Btilde,k)))
# }
# if(parallel==TRUE){
#   closeAllConnections()
#   cl <- makeCluster(detectCores())
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- 1:MaxLag
#   fe_call <- as.call( c(list (as.name("foreach"), i = i,.export=c("crossDist","ATilde")) ))
#   fe <- eval(fe_call)
#   Rstar <- fe %dopar% rstar(i)
#   stopCluster(cl)
#   quant <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],0.95))
#   pv <- sapply(1:MaxLag,FUN=function(j) mean(Rstar[[j]]>=quant[j]))
#   pvadj <- sapply(1:MaxLag,FUN=function(j) p.adjust(pv[j],method="fdr"))
#   cvadj <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],1-pvadj[j]))
#   res <- max(cvadj)
# }
# else {
#  Rstar <- lapply(1:MaxLag,FUN=function(j) rstar(j))
#  quant <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],0.95))
#  pv <- sapply(1:MaxLag,FUN=function(j) mean(Rstar[[j]]>=quant[j]))
#  pvadj <- sapply(1:MaxLag,FUN=function(j) p.adjust(pv[j],method="fdr"))
#  cvadj <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],1-pvadj[j]))
#  res <- max(cvadj)
# }
# return(res)
# }
