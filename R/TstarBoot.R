TstarBoot <- function(x,type,testType,p,b,parallel=FALSE){
 n <- length(x)
 if (is.matrix(x)) {
  if (!NCOL(x)==1) stop('Univariate time series only')
 } else {
  x <- c(x)
 }
 if(!all(is.finite(x))) stop('Missing or infitive values')
 if (!is.numeric(x))
     stop("'x' must be numeric")
 if((b==0) | missing(b)) stop('b must be grater than 0')
 A0 <- crossDist(x,lags=0)$A[[1]]
 Atilde0 <- ATilde(A0)
 MaxLag <- n-1
 test <- function(k){
  kern <- kernelFun(type,k/p)
  if (kern==0){
   d=rep(0,b)
  } else {
  cross <- crossDist(x,lags=k)
  A <- cross$A[[1]]
  B <- cross$B[[1]]
  Atilde <- ATilde(A)
  Btilde <- ATilde(B)
  bootCov = function(Atilde,Btilde,k){
   Wtstar <- rbind(rnorm(n-k))
    V <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))   
    return((n-k)*kern^2*V^2)
  }
  bootCor = function(Atilde,Btilde,k){
    Wtstar <- rbind(rnorm(n-k))
    dcov <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))
    dvarx <- sqrt(mean((Atilde0*Atilde0))*mean((Atilde0*Atilde0)))
    V <- dcov/sqrt(dvarx)
     return((n-k)*kern^2*V^2)
  }
  if (testType=="covariance"){
   d=replicate(b,bootCov(Atilde,Btilde,k))
  } else {
   d=replicate(b,bootCor(Atilde,Btilde,k))
  }
 }
d
}
if(parallel==TRUE){
  closeAllConnections()
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  clusterSetRNGStream(cl = cl, iseed = 9182)
  i <- seq_len(MaxLag)
  fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","crossDist","ATilde")) ))
  fe <- eval(fe_call)
  res <- fe %dopar% test(i)
  stopCluster(cl)
}
else {
 res <- rowSums(sapply(1:MaxLag,function(k) test(k)))
}
return(res)
}
