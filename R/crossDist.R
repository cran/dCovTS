crossDist <-
function(x,lags){
 if (is.matrix(x)==F) x <- as.matrix(x)
 n <- NROW(x)
 p <- NCOL(x)
 Y <- as.matrix(x[1:(n-lags),])
 X <- as.matrix(x[(1+lags):n,])
 A <- lapply(1:p,function(j) as.matrix(dist(X[,j],diag=TRUE,upper=TRUE)))
 B <- lapply(1:p,function(j) as.matrix(dist(Y[,j],diag=TRUE,upper=TRUE)))
 return(list(A=A,B=B))
}
