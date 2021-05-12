OrdinaryBoot2 <- function(x, type, p, b, parallel = FALSE) {

  if (is.vector(x) )  stop( 'Multivariate time series only.' )
  if ( !all( is.finite(x) ) )  stop( 'Missing or infitive values.' )
  if ( !is.numeric(x) )  stop( "'x' must be numeric." )

  n <- dim(x)[1]
  MaxLag <- n - 2

  test <- function(j) {
    kern <- kernelFun(type, j/p)
    if ( abs(kern) < 1e-16 ) {
     d <- numeric(b)
    } else {
      boot <- function(x, j) {
        ind <- sample(1:n, replace = TRUE)
        xStar <- x[ind, ]
        Rrm <- mADCF(xStar, j, unbiased = FALSE, output = FALSE)
        return( sum(Rrm^2) )
      }
      d <- (n - j) * kern^2 * replicate(b, boot(x, j) )
    }  ## if ( kern == 0 )
    d
  }

  if ( parallel ) {
    closeAllConnections()
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    clusterSetRNGStream(cl = cl, iseed = 9182)
    i <- seq_len(MaxLag)
    fe_call <- as.call( c(list (as.name("foreach"), i = i, .combine = "+",
                        .export = c("kernelFun", "mADCF"), .packages = "dcov" ) ) )
    fe <- eval(fe_call)
    res <- fe %dopar% test(i)
    stopCluster(cl)
  } else {
    res <- rowSums( sapply( 1:MaxLag, function(i) test(i) ) )
  }

  return(res)
}











# OrdinaryBoot2 <- function(x,type,p,b,parallel=FALSE){
#
#  if(is.vector(x))stop('Multivariate time series only.')
#  if(!all(is.finite(x))) stop('Missing or infitive values.')
#  if (!is.numeric(x)) stop("'x' must be numeric.")
#  n <- as.integer(NROW(x))
#  q <- as.integer(NCOL(x))
#  MaxLag <- n-2
#  test <- function(j){
#   kern <- kernelFun(type,j/p)
#   if (kern==0){
#    d=rep(0,b)
#   } else {
#   boot <- function(x,j){
#    ind <- sample(1:n,replace=T)
#    xStar <- x[ind,]
#    Rrm <- mADCF(xStar,j,unbiased=FALSE,output=FALSE)
#     return((n-j)*kern^2*sum(Rrm^2))
#   }
#   d <- replicate(b,boot(x,j))
#  }
#  d
# }
# if(parallel==TRUE){
#   closeAllConnections()
#   cl <- makeCluster(detectCores())
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- seq_len(MaxLag)
#   fe_call <- as.call(c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","mADCF","dcor")) ))
#   fe <- eval(fe_call)
#   res <- fe %dopar% test(i)
#   stopCluster(cl)
# }
# else {
#  res <- rowSums(sapply(1:MaxLag,function(i) test(i)))
# }
# return(res)
# }
