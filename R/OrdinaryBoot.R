OrdinaryBoot <- function(x, type, testType, p, b, parallel = FALSE) {

  n <- length(x)
  if ( is.matrix(x) ) {
    if ( dim(x)[2] > 1 ) {
     stop( "Univariate time series only." )
    } else  x <- x[, 1]
  }
  if ( !all(is.finite(x) ) )  stop( 'Missing or infitive values.' )
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  if( ( b == 0 )  |  missing(b) )  stop( 'b must be grater than 0.' )
  MaxLag <- n - 1

  if ( parallel ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- parallel::makePSOCKcluster( parallel::detectCores() )
    doParallel::registerDoParallel(cl)

    if ( testType == "covariance" ) {
      x_star <- replicate( b, Rfast2::Sample(x, n, replace = TRUE) )
      k <- 1:MaxLag
      d <- foreach(k = k, .combine = cbind, .packages = "dcov", .export = "kernelFun") %dopar% {
        kern <- kernelFun(type, k/p)
        if ( abs(kern) < 1e-16 ) {
          return( numeric(b) )
        } else {
          a <- numeric(b)
          xA <- x_star[ 1:(n - k), , drop = FALSE]
          xB <- x_star[ (1 + k):n, , drop = FALSE]
          for (j in 1:b)  a[j] <- dcov::dcov2d(xA[, j], xB[, j])
          return( (n - k) * kern^2 * a )
        } ## end if ( kern == 0 )
      } ## end foreach

    } else {
      x_star <- replicate( b, Rfast2::Sample(x, n, replace = TRUE) )
      k <- 1:MaxLag
      d <- foreach(k = k, .combine = cbind, .packages = "dcov", .export = "kernelFun") %dopar% {
        kern <- kernelFun(type, k/p)
        if ( abs(kern) < 1e-16 ) {
          return( numeric(b) )
        } else {
          a <- numeric(b)
          xA <- x_star[ 1:(n - k), , drop = FALSE]
          xB <- x_star[ (1 + k):n, , drop = FALSE]
          for (j in 1:b)  a[j] <- dcov::dcor2d(xA[, j], xB[, j])
          return( (n - k) * kern^2 * a )
        } ## end if ( kern == 0 )
      } ## end foreach

    } ## end if ( testType == "covariance" )
    parallel::stopCluster(cl)


  } else {
    if ( testType == "covariance" ) {
      x_star <- replicate( b, Rfast2::Sample(x, n, replace = TRUE) )
      d <- matrix(nrow = b, ncol = MaxLag)
      for (k in 1:MaxLag) {
        kern <- kernelFun(type, k/p)
        if ( abs(kern) < 1e-16 ) {
          d[, k] <- 0
        } else {
          a <- numeric(b)
          xA <- x_star[ 1:(n - k), , drop = FALSE]
          xB <- x_star[ (1 + k):n, , drop = FALSE]
          for (j in 1:b)  a[j] <- dcov::dcov2d(xA[, j], xB[, j])
          d[, k] <- (n - k) * kern^2 * a
        } ## end if ( kern == 0 )
      } ## end for (k in 1:MaxLag)

    } else {
      x_star <- replicate( b, Rfast2::Sample(x, n, replace = TRUE) )
      d <- matrix(nrow = b, ncol = MaxLag)
      for (k in 1:MaxLag) {
        kern <- kernelFun(type, k/p)
        if ( abs(kern) < 1e-16 ) {
          d[, k] <- 0
        } else {
          a <- numeric(b)
          xA <- x_star[ 1:(n - k), , drop = FALSE]
          xB <- x_star[ (1 + k):n, , drop = FALSE]
          for (j in 1:b)  a[j] <- dcov::dcor2d(xA[, j], xB[, j])
          d[, k] <- (n - k) * kern^2 * a
        } ## end if ( kern == 0 )
      } ## end for (k in 1:MaxLag)
    } ## end if ( tesType == "covariance" )

  } ## end if ( parallel )

  Rfast::rowsums(d)
}











# OrdinaryBoot <- function(x,type,testType,p,b,parallel=FALSE){
#  n <- length(x)
#  if (is.matrix(x)) {
#   if (!NCOL(x)==1) stop('Univariate time series only')
#  } else {
#   x <- c(x)
#  }
#  if(!all(is.finite(x))) stop('Missing or infitive values')
#  if (!is.numeric(x))
#      stop("'x' must be numeric")
#  if((b==0) | missing(b)) stop('b must be grater than 0')
#  MaxLag <- n-1
#  dcovFun <- function(x,k){
#      xA <- x[1:(n-k)]
#      xB <- x[(1+k):n]
#      return(dcov(xA,xB))
#  }
#  dcorFun <- function(x,k){
#      xA <- x[1:(n-k)]
#      xB <- x[(1+k):n]
#      return(dcor(xA,xB))
#  }
#  test <- function(k){
#   kern <- kernelFun(type,k/p)
#   if (kern==0){
#    d=rep(0,b)
#   } else {
#   if (testType=="covariance"){
#    d=replicate(b,expr={
#     x_star <- sample(x,replace=TRUE)
#     adcov <- dcovFun(x_star,k)
#     return((n-k)*kern^2*adcov^2)
#    })
#   } else {
#    d=replicate(b,expr={
#     x_star <- sample(x,replace=T)
#     adcor <- dcorFun(x_star,k)
#     return((n-k)*kern^2*adcor^2)
#   })
#   }
#  }
# d
# }
# if(parallel==TRUE){
#   closeAllConnections()
#   cl <- makeCluster(detectCores())
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- seq_len(MaxLag)
#   fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","dcor","dcov","bcdcor","ADCF","ADCV","dcovU")) ))
#   fe <- eval(fe_call)
#   res <- fe %dopar% test(i)
#   stopCluster(cl)
# }
# else {
#  res <- rowSums(sapply(1:MaxLag,function(k) test(k)))
# }
# return(res)
# }
