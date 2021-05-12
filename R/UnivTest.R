UnivTest <- function( x, type = c("truncated", "bartlett", "daniell", "QS", "parzen"),
                      testType = c("covariance","correlation"), p, b = 0, parallel = FALSE,
                      bootMethod = c("Wild Bootstrap", "Independent Bootstrap") ) {

  data.name <- deparse( substitute(x) )
  if ( is.matrix(x) ) {
    if ( dim(x)[2] > 1 ) {
     stop( "Univariate time series only." )
    } else  x <- x[, 1]
  }
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  if ( !all( is.finite(x) ) )  stop( "Missing or infitive values." )

  type <- match.arg(type)
  testType <- match.arg(testType)
  bootMethod <- match.arg(bootMethod)
  if ( missing(testType) )  method <- "covariance"
  if ( missing(bootMethod) )  method <- "Wild Bootstrap"

  n <- length(x)
  MaxLag <- n - 1

  adcv <- function(k, x, n) {
    xA <- x[1:(n - k)]
    xB <- x[(1 + k):n]
    adcv <- dcov::dcov2d(xA, xB)
    return(adcv)
  }
  adcf <- function(k, x, n) {
    xA <- x[1:(n - k)]
    xB <- x[(1 + k):n]
    adcf <- dcov::dcor2d(xA, xB)
    return(adcf)
  }

  ta <- numeric(MaxLag)
  if ( testType == "covariance" ) {
    for ( k in 1:MaxLag ) {
      kern <- kernelFun(type, k/p)
      if ( abs(kern) > 1e-16 )  ta[k] <- (n - k) * kern^2 * adcv(k, x, n)
    }
  } else {
    for ( k in 1:MaxLag ) {
      kern <- kernelFun(type, k/p)
      if ( abs(kern) > 1e-16 )  ta[k] <- (n - k) * kern^2 * adcf(k, x, n)
    }
  }

  method2 <- ifelse( ( testType == "covariance"), "Univariate test of independence based on distance covariance",
                                                  "Univariate test of independence based on distance correlation" )
  stat <- sum(ta)

  if ( !b == 0 ) {
    if ( bootMethod == "Wild Bootstrap" ) {
      Tnstar <- TstarBoot(x, type, testType, p, b, parallel)
    } else Tnstar <- OrdinaryBoot(x, type, testType, p, b, parallel)
    pvalue <- sum(Tnstar >= stat) / (b + 1)
  }

  p.value <- ifelse(b == 0, NA, pvalue)
  if ( b == 0 ) {
    Tnstar <- NULL
  } else  Tnstar <- Tnstar
  dataname <- paste( data.name, ",", " kernel type: ", type, ", bandwidth=", p,
                     ", replicates ", b, ", boot method: ", bootMethod, sep = "")

  names(stat) <- "Tn"
  e <- list(method = method2, statistic = stat, p.value = p.value, replicates = Tnstar,
            bootMethod = bootMethod, data.name = dataname)
  class(e) <- "htest"
  return(e)
}


















# UnivTest <- function (x, type = c("truncated", "bartlett", "daniell", "QS",
#  "parzen"), testType = c("covariance","correlation"), p, b = 0, parallel = FALSE, bootMethod=c("Wild Bootstrap","Independent Bootstrap"))
# {
#     type <- match.arg(type)
#     testType <- match.arg(testType)
#     bootMethod <- match.arg(bootMethod)
#     if (missing(testType)) method="covariance"
#     if (missing(bootMethod)) method="Wild Bootstrap"
#     data.name <- deparse(substitute(x))
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
#
#     MaxLag <- n - 1
#     adcv <- function(k,x){
#         n <- length(x)
#         xA <- x[1:(n-k)]
#         xB <- x[(1+k):n]
#         adcv <- dcov(xA,xB)
#       return(adcv)
#      }
#     adcf <- function(k,x){
#         n <- length(x)
#         xA <- x[1:(n-k)]
#         xB <- x[(1+k):n]
#         adcf <- dcor(xA,xB)
#     return(adcf)
#     }
#     t = rep(0, MaxLag)
#     for (k in 1:MaxLag) {
#         kern <- kernelFun(type, k/p)
#         if (kern != 0) {
#                 if(testType=="covariance"){
#                   t[k] <- (n - k) * kern^2 * adcv(k,x)^2
#                 } else {
#                   t[k] <- (n - k) * kern^2 * adcf(k,x)^2
#                 }
#         }
#     }
#     method2 = ifelse((testType=="covariance"),"Univariate test of independence based on distance covariance",
# "Univariate test of independence based on distance correlation")
#     stat <- sum(t)
#     if (!b == 0) {
#         if (bootMethod=="Wild Bootstrap"){
#          Tnstar <- TstarBoot(x, type, testType, p, b, parallel)
#         }
#         else {
#          Tnstar <- OrdinaryBoot(x, type, testType, p, b, parallel)
#         }
#         pvalue <- sum(Tnstar >= stat)/(b + 1)
#     }
#     p.value <- ifelse(b == 0, NA, pvalue)
#     if (b == 0) {
#         Tnstar <- NULL
#     }
#     else Tnstar <- Tnstar
#     dataname <- paste(data.name, ",", " kernel type: ", type,
#         ", bandwidth=", p, ", replicates ", b, ", boot method: ", bootMethod, sep = "")
#
#     names(stat) <- "Tn"
#     e = list(method = method2, statistic = stat, p.value = p.value, replicates = Tnstar, bootMethod = bootMethod,
#         data.name = dataname)
#     class(e) = "htest"
#     return(e)
# }


