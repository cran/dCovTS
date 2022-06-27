OrdinaryBootCV <- function(n, MaxLag, alpha = 0.05, b, parallel = FALSE) {

   if (n < 500) {
     x <- rnorm(n)
   } else x <- Rfast::Rnorm(n)
   xStar <- replicate( b, Rfast2::Sample(x, n, replace = TRUE) )
   Rstar <- matrix(nrow = b, ncol = MaxLag )

   if ( parallel ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
     requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
     closeAllConnections()
     cl <- parallel::makePSOCKcluster( parallel::detectCores() )
     doParallel::registerDoParallel(cl)
     rstarb <- numeric(b)
     k <- 1:MaxLag
     Rstar <- foreach(k = k, .combine = cbind, .packages = "dcov" ) %dopar% {
       xA <- xStar[1:(n - k), ]
       xB <- xStar[(1 + k):n, ]
       for (i in 1:b)  rstarb[i] <- dcov::dcor2d(xA[, i], xB[, i])
       return( rstarb )
     }
     parallel::stopCluster(cl)

   } else {
     for ( k in 1:MaxLag ) {
       xA <- xStar[1:(n - k), ]
       xB <- xStar[(1 + k):n, ]
       for (i in 1:b)  Rstar[i, k] <- dcov::dcor2d(xA[, i], xB[, i])
     }
   }

   Rstar <- sqrt(Rstar)
   quant <- Rfast2::colQuantile(Rstar, probs = 1 - alpha)
   pv <- pvadj <- Rfast::colmeans( Rfast::eachrow( Rstar, quant, oper = "-") >= 0)
   for (i in 1:k)  pvadj[i] <- p.adjust(pv[i], method = "fdr")
   cvadj <- Rfast2::colQuantile(Rstar, 1 - pvadj)
   max(cvadj)
}












# OrdinaryBootCV <- function(n,MaxLag,b,parallel=FALSE){
#  x <- rnorm(n)
#  if (missing(MaxLag) || MaxLag < 0)
#       stop("'MaxLag' must be greater than 1")
#  dcorFun <- function(x,k){
#      xA <- x[1:(n-k)]
#      xB <- x[(1+k):n]
#      return(dcor(xA,xB))
#  }
#  rstar <- function(k){
#   Rstark <- function(k){
#    xStar <- sample(x,replace=TRUE)
#    dcor <- dcorFun(xStar,k)
#    return(dcor)
#   }
#  return(replicate(b,Rstark(k)))
# }
# if(parallel==TRUE){
#   closeAllConnections()
# #  cl <- makeCluster(detectCores())
#   cl <- makeCluster(2)
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- 1:MaxLag
#   fe_call <- as.call( c(list (as.name("foreach"), i = i,.export=c("dcor","ADCF")) ))
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
