OrdinaryBoot <- function(x,type,testType,p,b,parallel=FALSE){
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
 MaxLag <- n-1
 dcovFun <- function(x,k){
     xA <- x[1:(n-k)]
     xB <- x[(1+k):n]
     return(dcov(xA,xB)) 
 }
 dcorFun <- function(x,k){
     xA <- x[1:(n-k)]
     xB <- x[(1+k):n]
     return(dcor(xA,xB))
 } 
 test <- function(k){
  kern <- kernelFun(type,k/p)
  if (kern==0){
   d=rep(0,b)
  } else {
  if (testType=="covariance"){
   d=replicate(b,expr={
    x_star <- sample(x,replace=TRUE)
    adcov <- dcovFun(x_star,k)
    return((n-k)*kern^2*adcov^2)
   })
  } else {
   d=replicate(b,expr={
    x_star <- sample(x,replace=T)
    adcor <- dcorFun(x_star,k)
    return((n-k)*kern^2*adcor^2)
  })
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
  fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","dcor","dcov","bcdcor","ADCF","ADCV","dcovU")) ))
  fe <- eval(fe_call)
  res <- fe %dopar% test(i)
  stopCluster(cl)
}
else {
 res <- rowSums(sapply(1:MaxLag,function(k) test(k)))
}
return(res)
}

