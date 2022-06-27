TstarBoot <- function(x, type, testType, p, b, parallel = FALSE) {

  if ( is.matrix(x) ) {
    if ( dim(x)[2] > 1 ) {
      stop( 'Univariate time series only.' )
    } else  x <- x[, 1]
  }
  if ( !all( is.finite(x) ) )  stop( 'Missing or infitive values.' )
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  if ( (b == 0)  |  missing(b) )  stop( 'b must be grater than 0.' )

  n <- length(x)
  MaxLag <- n - 1
  Atilde0 <- ATilde( Rfast::vecdist(x) )
  dvarx <- mean(Atilde0 * Atilde0)

  bootCov <- function(Atilde, Btilde, k, b) {
    if ( b %% 2 == 0 ) {
      b1 <- b2 <- b/2
    } else {
      b1 <- floor(b/2)
      b2 <- ceiling(b/2)
    }
    com <- Atilde * Btilde
    Wtstar <- Rfast::matrnorm(b1, n - k) ## matrix of Z variables with b/2 rows and n - k columns
    V <- Rfast::rowsums( Wtstar %*% com * Wtstar ) / (n - k)  ## this is in fact a Mahalanobis distance
    Wtstar <- Rfast::matrnorm(b2, n - k) #
    V2 <- Rfast::rowsums( Wtstar %*% com * Wtstar ) / (n - k)
    V <- c(V, V2)
    return( kern^2 * V )
  }

  bootCor <- function(Atilde, Btilde, k, b, dvarx) {
    if ( b %% 2 == 0 ) {
      b1 <- b2 <- b/2
    } else {
      b1 <- floor(b/2)
      b2 <- ceiling(b/2)
    }
    com <- Atilde * Btilde
    Wtstar <- Rfast::matrnorm(b1, n - k) ## matrix of Z variables with b/2 rows and n - k columns
    V <- Rfast::rowsums( Wtstar %*% com * Wtstar ) / dvarx  ## this is in fact a Mahalanobis distance
    Wtstar <- Rfast::matrnorm(b2, n - k)
    V2 <- Rfast::rowsums( Wtstar %*% com * Wtstar ) / dvarx
    V <- c(V, V2)
    return( kern^2 * V / (n - k) )
  }

  if ( parallel ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- parallel::makePSOCKcluster( parallel::detectCores() )
    doParallel::registerDoParallel(cl)

    k <- 1:MaxLag
    d <- foreach(k = k, .combine = cbind, .export = c("kernelFun", "crossDist", "ATilde", "bootCov", "bootCor"), .packages = "Rfast") %dopar% {
      kern <- kernelFun(type, k/p)
      if ( abs(kern) < 1e-16 ) {
        return( numeric(b) )
      } else {
        cross <- crossDist(x, lags = k)
        Atilde <- ATilde( cross$A[[ 1 ]] )
        Btilde <- ATilde( cross$B[[ 1 ]] )
        if ( testType == "covariance" ) {
          return( bootCov(Atilde, Btilde, k, b) )
        } else  return( bootCor(Atilde, Btilde, k, b, dvarx) )
      } ##  end  if ( kern == 0 )
    }  ## end foreach
    parallel::stopCluster(cl)

  } else {

    d <- matrix( nrow = b, ncol = MaxLag)
    for ( k in 1:MaxLag ) {
      kern <- kernelFun(type, k/p)
      if ( abs(kern) < 1e-16 ) {
        d[, k] <- numeric(b)
      } else {
        cross <- crossDist(x, lags = k)
        Atilde <- ATilde( cross$A[[ 1 ]] )
        Btilde <- ATilde( cross$B[[ 1 ]] )
        if ( testType == "covariance" ) {
          d[, k] <- bootCov(Atilde, Btilde, k, b)
        } else  d[, k] <-  bootCor(Atilde, Btilde, k, b, dvarx)
      }  ##  end  if ( kern == 0 )
    }  ## end  for ( k in 1:MaxLag )

  } ##  end if ( parallel )

  Rfast::rowsums(d)
}










# TstarBoot <- function(x,type,testType,p,b,parallel=FALSE){
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
#  A0 <- crossDist(x,lags=0)$A[[1]]
#  Atilde0 <- ATilde(A0)
#  MaxLag <- n-1
#  test <- function(k){
#   kern <- kernelFun(type,k/p)
#   if (kern==0){
#    d=rep(0,b)
#   } else {
#   cross <- crossDist(x,lags=k)
#   A <- cross$A[[1]]
#   B <- cross$B[[1]]
#   Atilde <- ATilde(A)
#   Btilde <- ATilde(B)
#   bootCov = function(Atilde,Btilde,k){
#    Wtstar <- rbind(rnorm(n-k))
#     V <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))
#     return((n-k)*kern^2*V^2)
#   }
#   bootCor = function(Atilde,Btilde,k){
#     Wtstar <- rbind(rnorm(n-k))
#     dcov <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))
#     dvarx <- sqrt(mean((Atilde0*Atilde0))*mean((Atilde0*Atilde0)))
#     V <- dcov/sqrt(dvarx)
#      return((n-k)*kern^2*V^2)
#   }
#   if (testType=="covariance"){
#    d=replicate(b,bootCov(Atilde,Btilde,k))
#   } else {
#    d=replicate(b,bootCor(Atilde,Btilde,k))
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
#   fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","crossDist","ATilde")) ))
#   fe <- eval(fe_call)
#   res <- fe %dopar% test(i)
#   stopCluster(cl)
# }
# else {
#  res <- rowSums(sapply(1:MaxLag,function(k) test(k)))
# }
# return(res)
# }
