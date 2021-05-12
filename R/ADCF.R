ADCF <- function (x, MaxLag = 15, unbiased = FALSE) {

  n <- length(x)
  if ( is.matrix(x) ) {
    if ( dim(x)[2] > 1 ) {
     stop( "Univariate time series only." )
    } else  x <- x[, 1]
  }
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  if ( !all( is.finite(x) ) )  stop( "Missing or infitive values." )
  if ( ( MaxLag < 0 )  |  ( MaxLag > n - 1 ) )  stop( "'MaxLag' must be in the range of 0 and (n-1)." )
  if ( ( unbiased )  &&  ( n < 3 ) )  stop( "Sample size must be larger than 3." )

  adcf <- matrix(0, 1, MaxLag + 1, dimnames = list( "  ", paste( "lag", 0:MaxLag) ) )

  if ( unbiased ) {
    for ( k in 0:MaxLag ) {
      xA <- x[1:(n - k)]
      xB <- x[(1 + k):n]
      adcf[, (k + 1)] <- dcov::dcor2d(xA, xB, type = "U")
    }

  } else {
    for ( k in 0:MaxLag ) {
      xA <- x[1:(n - k)]
      xB <- x[(1 + k):n]
      adcf[, (k + 1)] <- sqrt( dcov::dcor2d(xA, xB) )
    }
  }

  return(adcf)
}










# ADCF <- function (x, MaxLag=15,unbiased=FALSE) {
#     n <- length(x)
#     if (is.matrix(x)) {
#         if (!NCOL(x) == 1)
#             stop("Univariate time series only")
#     }
#     else {
#         x <- c(x)
#     }
#     if (!is.numeric(x))
#         stop("'x' must be numeric")
#     if (!all(is.finite(x)))
#         stop("Missing or infitive values")
#     if((MaxLag<0) | (MaxLag>(n-1)))stop("'MaxLag' must be in the range of 0 and (n-1)")
#     adcf <- matrix(0, 1,MaxLag+1,dimnames=list("  ",paste("lag",0:MaxLag)))
#     for (k in seq(0,MaxLag,by=1)){
#      xA <- x[1:(n-k)]
#      xB <- x[(1+k):n]
#      if (unbiased) {
#       adcf[,(k+1)] <- round(bcdcor(xA,xB),4)
#      }
#      else {
#       adcf[,(k+1)] <- round(dcor(xA,xB),4)
#      }
#     }
#     return(adcf)
# }






