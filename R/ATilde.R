ATilde <- function(A) {

  n <- dim(A)[1]
  total.mean <- sum(A)/n^2
  m1 <- Rfast::rowmeans(A)
  m2 <- Rfast::colmeans(A)
  A2 <- t( t( t(A - m1) - m2 ) ) + total.mean
  A2
}








# ATilde <- function(A){
#  n <- NROW(A)
#  total.mean <- sum(A)/(n^2)
#  m1 <- sapply(1:n, FUN=function(j) mean(A[j,]))
#  m2 <- sapply(1:n, FUN=function(j) mean(A[,j]))
#  A2 <- sapply(1:n,1:n,FUN=function(i,j) (A[i,j]-m1[i]-m2[j]+total.mean))
#  return(A2)
# }
