mADCV <- function(x, lags, unbiased = FALSE, output = TRUE) {

  if ( !is.matrix(x) )  x <- as.matrix(x)
  if ( (is.data.frame(x) )  |  ( is.matrix(x) ) ) {
    if ( dim(x)[2] == 1 )  stop( 'Only multivariate time series with dimension d>2.' )
    if ( !is.ts(x) )  x <- as.ts(x)
  }
  if ( !is.numeric(x) )  stop( "'x' must be numeric" )
  if ( !all( is.finite(x) ) )  stop( 'Missing or infitive values.' )
  if ( missing(lags) )  stop( "'lags' is missing with no default." )

  dm <- dim(x)
  n <- dm[1]   ;   d <- dm[2]

  if ( length(lags) == 1 ) {

    if ( ( lags >= 0 )  &&  ( lags <= (n - 1) ) ) {
      X <- x[(1 + lags):n, ]
      Y <- x[1:(n - lags), ]
    } else if ( ( lags >= (-n + 1) )  &&  ( lags < 0 ) ) {
      X <- x[1:(n + lags), ]
      Y <- x[(1 - lags):n, ]
    } else  stop( "lags must be in the range of -(n-1) and (n-1)." )

    cross.adcv <- matrix(NA, d, d)
    if ( unbiased ) {
      for ( i in 1:d )  cross.adcv[i, ] <- dcov::mdcov(X[, i], Y, type = "U")
    } else {
      for ( i in 1:d )  cross.adcv[i, ] <- dcov::mdcov(X[, i], Y)
      cross.adcv <- sqrt(cross.adcv)
    }

    if ( output ) {
      cat( "Distance Covariance Matrix at lag: ", lags, "\n" )
      print(cross.adcv)
    }

  } else if ( length(lags) > 1 ) {
    len <- length(lags)
    cross.adcv <- array( dim = c(d, d, len))

    for ( j in 1:len ) {
      cross <- matrix(NA, d, d)
      lag <- lags[j]

      if ( ( lag >= 0 )  &&  ( lag <= (n - 1) ) ) {
        X <- x[(1 + lag):n, ]
        Y <- x[1:(n - lag), ]
      } else if ( ( lag >= (-n + 1) )  &&  ( lag < 0 ) ) {
        X <- x[1:(n + lag), ]
        Y <- x[(1 - lag):n, ]
      } else  stop( "lags must be in the range of -(n-1) and (n-1)." )

      if ( unbiased ) {
        for ( i in 1:d )  cross[i, ] <- dcov::mdcov(X[, i], Y, type = "U")
      } else {
        for ( i in 1:d )  cross[i, ] <- dcov::mdcov(X[, i], Y)
        cross <- sqrt(cross)
      }
      cross.adcv[, , j] <- cross
    } ## end  for ( i in 1:len ) 

    dimnames(cross.adcv)[[ 3 ]] <- paste( "Distance Covariance Matrix at lag: ", lags) 
    if ( output )  print(cross.adcv)

  } ## end  if ( length(lags) == 1 ) 


  mADCV <- cross.adcv
}











# mADCV <- function(x,lags,unbiased=FALSE,output=TRUE){
#  if(!is.matrix(x)) x <- as.matrix(x)
#  if ((is.data.frame(x)) | (is.matrix(x))){
#   if(NCOL(x)==1)stop('Only multivariate time series with dimension d>2')
#   if(!is.ts(x)) x <- as.ts(x)
#  }
#  if (!is.numeric(x))
#      stop("'x' must be numeric")
#  if(!all(is.finite(x))) stop('Missing or infitive values')
#  if (missing(lags))
#        stop("'lags' is missing with no default")
#  n <- as.integer(NROW(x))
#  d <- as.integer(NCOL(x))
#  if ((lags >= 0) && (lags <= (n-1))){
#   X <- rbind(x[(1+lags):n,])
#   Y <- rbind(x[1:(n-lags),])
#  } else if ((lags >= (-n+1)) && (lags <0)){
#   X <- rbind(x[1:(n+lags),])
#   Y <- rbind(x[(1-lags):n,])
#  } else {
#   stop("'lags' must be in the range of -(n-1) and (n-1)")
#  }
#  cross.adcv <- matrix(NA,d,d)
#  for (i in 1:d){
#   for (j in 1:d){
#     if (unbiased){
#   cross.adcv[i,j] <- dcovU(X[,i],Y[,j])
#  }
#  else {
#    cross.adcv[i,j] <- dcov(X[,i],Y[,j])
#  }
#   }
#  }
#  if(output){
#  cat("Distance Covariance Matrix at lag: ", lags, "\n")
#  print(cross.adcv)
#  }
#  mADCV <- cross.adcv
# }
