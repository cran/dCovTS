crossDist <- function(x, lags) {

  if ( !is.matrix(x) )  x <- as.matrix(x)
  dm <- dim(x)
  n <- dm[1]   ;   p <- dm[2]
  Y <- x[1:(n - lags), , drop = FALSE]
  X <- x[(1 + lags):n, , drop = FALSE]
  a <- paste("A", 1:p)
  A <- sapply(a, function(x) NULL)
  b <- paste("B", 1:p)
  B <- sapply(b, function(x) NULL)

  if ( lags > 0 ) {
    for (j in 1:p) {
      A[[ j ]] <- Rfast::vecdist( X[, j] )
      B[[ j ]] <- Rfast::vecdist( Y[, j] )
    }
  } else {
    for (j in 1:p)  A[[ j ]] <- B[[ j]] <- Rfast::vecdist( X[, j] )
  }

  list(A = A, B = B)
}






# crossDist <- function(x,lags){
#  if (is.matrix(x)==F) x <- as.matrix(x)
#  n <- NROW(x)
#  p <- NCOL(x)
#  Y <- as.matrix(x[1:(n-lags),])
#  X <- as.matrix(x[(1+lags):n,])
#  A <- lapply(1:p,function(j) as.matrix(dist(X[,j],diag=TRUE,upper=TRUE)))
#  B <- lapply(1:p,function(j) as.matrix(dist(Y[,j],diag=TRUE,upper=TRUE)))
#  return(list(A=A,B=B))
# }
