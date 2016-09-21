OrdinaryBootCV <- function(n,MaxLag,b,parallel=FALSE){
 x <- rnorm(n)
 if (missing(MaxLag) || MaxLag < 0)
      stop("'MaxLag' must be greater than 1")
 dcorFun <- function(x,k){
     xA <- x[1:(n-k)]
     xB <- x[(1+k):n]
     return(dcor(xA,xB))
 }
 rstar <- function(k){
  Rstark <- function(k){
   xStar <- sample(x,replace=TRUE)
   dcor <- dcorFun(xStar,k)
   return(dcor)
  }
 return(replicate(b,Rstark(k)))
}
if(parallel==TRUE){
  closeAllConnections()
#  cl <- makeCluster(detectCores())
  cl <- makeCluster(2)
  registerDoParallel(cl)
  clusterSetRNGStream(cl = cl, iseed = 9182)
  i <- 1:MaxLag
  fe_call <- as.call( c(list (as.name("foreach"), i = i,.export=c("dcor","ADCF")) ))
  fe <- eval(fe_call)
  Rstar <- fe %dopar% rstar(i)
  stopCluster(cl)
  quant <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],0.95))
  pv <- sapply(1:MaxLag,FUN=function(j) mean(Rstar[[j]]>=quant[j]))
  pvadj <- sapply(1:MaxLag,FUN=function(j) p.adjust(pv[j],method="fdr"))
  cvadj <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],1-pvadj[j]))
  res <- max(cvadj)
}
else {
 Rstar <- lapply(1:MaxLag,FUN=function(j) rstar(j))
 quant <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],0.95))
 pv <- sapply(1:MaxLag,FUN=function(j) mean(Rstar[[j]]>=quant[j]))
 pvadj <- sapply(1:MaxLag,FUN=function(j) p.adjust(pv[j],method="fdr"))
 cvadj <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],1-pvadj[j]))
 res <- max(cvadj)
}
return(res)
}

