mADCVtest <- function(x, type = c("truncated", "bartlett", "daniell", "QS", "parzen"), p, b = 0,
                      parallel = FALSE, bootMethod = c( "Wild Bootstrap","Independent Bootstrap" ) ) {

  type <- match.arg(type)
  data.name <- deparse( substitute(x) )

  if ( !is.matrix(x) )  stop( 'Only multivariate time series with dimension d>2.' )
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  if( !all( is.finite(x) ) )  stop( 'Missing or infitive values.' )
  n <- dim(x)[1]

  bootMethod <- match.arg(bootMethod)
  if ( missing(bootMethod) )  method <- "Wild Bootstrap"
  MaxLag <- n - 2
  ta <- numeric(MaxLag)

  for ( k in 1:MaxLag ) {
    kern <- kernelFun(type, k/p)
    if ( abs(kern) > 1e-16 )  ta[k] <- (n - k) * kern^2 * sum( mADCV(x, lags = k, unbiased = FALSE, output = FALSE)^2 )
  }
  stat <- sum(ta)

  if ( !b == 0 ) {
    if ( bootMethod == "Wild Bootstrap" ) {
      Tnstar <- TstarBoot1(x, type, p, b, parallel)
    } else  Tnstar <- OrdinaryBoot1(x, type, p, b, parallel)
      pvalue <- sum(Tnstar >= stat) / (b + 1)
  }
  p.value <- ifelse(b == 0, NA, pvalue)

  if ( b == 0 ) {
    Tnstar <- NULL
  } else Tnstar <- Tnstar
  dataname <- paste( data.name,","," kernel type: ", type,", bandwidth=",p, ", replicates ", b,
                     ", boot method: ", bootMethod, sep = "")
  names(stat) <- "Tn"
  e <- list(method = paste( "Multivariate test of independence based on distance covariance", sep = "" ),
            statistic = stat, p.value = p.value, replicates = Tnstar, bootMethod = bootMethod, data.name = dataname )
  class(e) <- "htest"

  return(e)
}


















# mADCVtest <- function(x,type = c("truncated", "bartlett", "daniell", "QS",
#  "parzen"), p, b = 0, parallel = FALSE, bootMethod=c("Wild Bootstrap","Independent Bootstrap"))
# {
#  type <- match.arg(type)
#  data.name <- deparse(substitute(x))
#  if (!is.matrix(x))stop('Only multivariate time series with dimension d>2')
#  if (!is.numeric(x))
#      stop("'x' must be numeric")
#  if(!all(is.finite(x))) stop('Missing or infitive values')
#  bootMethod <- match.arg(bootMethod)
#  if (missing(bootMethod)) method="Wild Bootstrap"
#  n <- as.integer(NROW(x))
#  q <- as.integer(NCOL(x))
#  MaxLag <- n-1
#  t <- rep(0,MaxLag)
#  for(k in 1:MaxLag){
#   kern <- kernelFun(type,k/p)
#   if (kern !=0){
#     t[k] <- (n-k)*kern^2*sum(mADCV(x,lags=k,unbiased=FALSE,output=FALSE)^2)
#   }
#  }
#  stat <- sum(t)
#  if(!b==0){
#   if (bootMethod=="Wild Bootstrap"){
#    Tnstar <- TstarBoot1(x,type,p,b,parallel)
#   }
#   else {
#    Tnstar <- OrdinaryBoot1(x,type,p,b,parallel)
#   }
#  pvalue <- sum(Tnstar>=stat)/(b+1)
#  }
#  p.value <- ifelse(b == 0,NA,pvalue)
#   if(b==0){
#    Tnstar <- NULL
#   } else Tnstar <- Tnstar
#  dataname <- paste(data.name,","," kernel type: ", type,", bandwidth=",p, ", replicates ", b,  ", boot method: ", bootMethod, sep = "")
#  names(stat) <- "Tn"
#    e=list(method = paste("Multivariate test of independence based on distance covariance", sep = ""),
#         statistic = stat, p.value = p.value, replicates=Tnstar,bootMethod=bootMethod,data.name=dataname)
#  class(e) <- "htest"
#  return(e)
# }
