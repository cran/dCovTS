TstarBoot2 <- function(x, type, p, b, parallel = FALSE) {

  if ( is.vector(x) )  stop( 'Multivariate time series only.' )
  if ( !all( is.finite(x) ) )  stop( 'Missing or infitive values.' )
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  dm <- dim(x)
  n <- dm[1]   ;   qa <- dm[2]
  A0 <- crossDist(x, 0)$A
  Atilde0 <- list()
  for ( i in 1:qa )  Atilde0[[ i ]] <- ATilde(A0[[ i ]] )
  
  boot <- function(Atilde0, Atilde, Btilde, kern, j, b) {
    ela <- expand.grid(1:qa, 1:qa)
    dvarx <- numeric( dim(ela)[1] )
    for ( i in 1:dim(ela)[1] ) {
      dvarx[i] <- mean( Atilde0[[ ela[i, 1] ]] ^2 ) * mean( Atilde0[[ ela[i, 2] ]]^2 )
    }
    dvarx <- sqrt(dvarx)
    if ( b %% 2 == 0 ) {
      b1 <- b2 <- b/2
    } else {
      b1 <- floor(b/2)
      b2 <- ceiling(b/2)
    }
    Rrm <- matrix(nrow = b1, ncol = dim(ela)[1])
    Wtstar <- Rfast::matrnorm(b1, n - j) ## matrix of Z variables with b/2 rows and n - k columns
    for ( i in 1:dim(ela)[1] ) {
      com <- Atilde[[ ela[i, 1] ]] * Btilde[[ ela[i, 2] ]]
      dcova <- Rfast::rowsums( Wtstar %*% com * Wtstar )
      Rrm[, i] <- dcova/dvarx[i]
    }
    a1 <- kern^2 * Rfast::rowsums(Rrm) / (n - j)
    Rrm <- matrix(nrow = b2, ncol = dim(ela)[1])   
    Wtstar <- Rfast::matrnorm(b2, n - j) ## matrix of Z variables with b/2 rows and n - k columns
    for ( i in 1:dim(ela)[1] ) {
      com <- Atilde[[ ela[i, 1] ]] * Btilde[[ ela[i, 2] ]]
      dcova <- Rfast::rowsums( Wtstar %*% com * Wtstar )
      Rrm[, i] <- dcova/dvarx[i]
    }
    a2 <- kern^2 * Rfast::rowsums(Rrm) / (n - j)
    return( c(a1, a2) )
  }  ##  end boot() function

  MaxLag <- n - 2
  test <- function(j) {
    kern <- kernelFun(type, j/p)
    if ( abs(kern) < 1e-16 ) {
      d <- numeric(b)
    } else {
      A <- crossDist(x, j)$A
      B <- crossDist(x, j)$B
      Atilde <- Btilde <- list()
      for ( i in 1:qa ) {
        Atilde[[ i ]] <- ATilde(A[[ i ]] )
        Btilde[[ i ]] <- ATilde(B[[ i ]] )
      }
      d <- boot(Atilde0, Atilde, Btilde, kern, j, b)
    }  ##  end if (kern == 0)
    d
  } ##  end test() function

  if ( parallel ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- makeCluster( detectCores() )
    registerDoParallel(cl)
    clusterSetRNGStream(cl = cl, iseed = 9182)
    i <- seq_len(MaxLag)
    fe_call <- as.call( c(list( as.name("foreach"), i = i, .combine="+",
               .export = c("kernelFun", "crossDist", "ATilde", "boot"), .packages = "Rfast" ) ) )
    fe <- eval(fe_call)
    res <- fe %dopar% test(i)
    stopCluster(cl)
  } else {
    res <- rowSums( sapply( 1:MaxLag, function(i) test(i) ) )
  }

  return(res)
}










# TstarBoot2 <- function(x,type,p,b,parallel=FALSE){
#  if(is.vector(x))stop('Multivariate time series only')
#  if(!all(is.finite(x))) stop('Missing or infitive values')
#  if (!is.numeric(x)) stop("'x' must be numeric")
#  n <- as.integer(NROW(x))
#  q <- as.integer(NCOL(x))
#  A0 <- crossDist(x,0)$A
#  Atilde0 <- lapply(1:q, FUN=function(j) ATilde(A0[[j]]))
#  MaxLag <- n-2
#  test <- function(j){
#   kern <- kernelFun(type,j/p)
#   if (kern==0){
#    d=rep(0,b)
#   } else {
#   A <- crossDist(x,j)$A
#   B <- crossDist(x,j)$B
#   Atilde <- lapply(1:q, FUN=function(i) ATilde(A[[i]]))
#   Btilde <- lapply(1:q, FUN=function(i) ATilde(B[[i]]))
#   boot <- function(Atilde0,Atilde,Btilde,j){
#    Wtstar <- rbind(rnorm(n-j))
#    Rrm <- matrix(NA,q,q)
#     for (l in 1:q){
#      for (f in 1:q){
#         dcov <- sqrt((Wtstar%*%(Atilde[[l]]*Btilde[[f]])%*%t(Wtstar))/((n-j)^2))
#         dvarx <- sqrt(mean((Atilde0[[l]]*Atilde0[[l]]))*mean((Atilde0[[f]]*Atilde0[[f]])))
#         Rrm[l,f] <- dcov/sqrt(dvarx)
#      }
#    }
#     return((n-j)*kern^2*sum(Rrm^2))
#   }
#   d <- replicate(b,boot(Atilde0,Atilde,Btilde,j))
#  }
#  d
# }
# if(parallel==TRUE){
#   closeAllConnections()
#   cl <- makeCluster(detectCores())
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- seq_len(MaxLag)
#   fe_call <- as.call(c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","crossDist","ATilde")) ))
#   fe <- eval(fe_call)
#   res <- fe %dopar% test(i)
#   stopCluster(cl)
# }
# else {
#  res <- rowSums(sapply(1:MaxLag,function(i) test(i)))
# }
# return(res)
# }


