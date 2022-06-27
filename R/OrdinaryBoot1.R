OrdinaryBoot1 <- function(x, type, p, b, parallel = FALSE) {

  if ( is.vector(x) )  stop( 'Multivariate time series only.' )
  if ( !all( is.finite(x) ) )  stop( 'Missing or infitive values.' )
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )
  n <- dim(x)[1]
  MaxLag <- n - 2

  boot <- function(x, j) {
    ind <- Rfast2::Sample.int(n, n, replace = TRUE )
    xStar <- x[ind, ]
    Vrm <- mADCV(xStar, j, unbiased = FALSE, output = FALSE)
    sum(Vrm^2)
  }

  test <- function(j) {
    kern <- kernelFun(type, j/p)
    if ( abs(kern) < 1e-16 ) {
      d <- numeric(b)
    } else  d <- (n - j) * kern^2 * replicate(b, boot(x, j) )
    d
  }

  if ( parallel ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    clusterSetRNGStream(cl = cl, iseed = 9182)
    i <- 1:MaxLag
    fe_call <- as.call( c( list( as.name("foreach"), i = i, .combine = "+",
               .export = c("kernelFun", "mADCV", "boot", "Sample.int"), .packages = c("Rfast2", "dcov") ) ) )
    fe <- eval(fe_call)
    res <- fe %dopar% test(i)
    stopCluster(cl)

  } else  res <- rowSums( sapply( 1:MaxLag, function(i) test(i) ) )

  res
}












# OrdinaryBoot1 <- function(x,type,p,b,parallel=FALSE){
#  if(is.vector(x))stop('Multivariate time series only')
#  if(!all(is.finite(x))) stop('Missing or infitive values')
#  if (!is.numeric(x)) stop("'x' must be numeric")
#  n <- as.integer(NROW(x))
#  q <- as.integer(NCOL(x))
#  MaxLag <- n-2
#  test <- function(j){
#   kern <- kernelFun(type,j/p)
#   if (kern==0){
#    d=rep(0,b)
#   } else {
#   boot = function(x,j){
#    ind <- sample(1:n,replace=T)
#    xStar <- x[ind,]
#    Vrm <- mADCV(xStar,j,unbiased=FALSE,output=FALSE)
#    res <- (n-j)*kern^2*sum(Vrm^2)
#    return(res)
#   }
#   d=replicate(b,boot(x,j))
#  }
#  d
# }
# if(parallel==TRUE){
#   closeAllConnections()
#   cl <- makeCluster(detectCores())
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- seq_len(MaxLag)
#   fe_call <- as.call(c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","mADCV","dcov")) ))
#   fe <- eval(fe_call)
#   res <- fe %dopar% test(i)
#   stopCluster(cl)
# }
# else {
#  res <- rowSums(sapply(1:MaxLag,function(i) test(i)))
# }
# return(res)
# }
