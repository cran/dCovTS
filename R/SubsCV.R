SubsCV <- function(x, MaxLag, alpha, parallel = FALSE) {

  n <- length(x)
  if ( ( n - MaxLag < 0 ) || ( n - MaxLag < 4 ) || ( n - MaxLag <= 25 ) )  stop( "Give bigger sample size n." )

  dcorFun <- function(x,k,b){
    xA <- x[1:(b - k)]
    xB <- x[(1 + k):b]
    sqrt( dcov::dcor2d(xA, xB) )
  }

  test <- function(x, j, b) {
    qa <- n - b + 1
    qa <- ifelse( qa < 0, abs(qa), qa)
    Xb <- sapply( 1:qa, function(i) x[i:(i + b - 1)] )
    Vstar <- unlist( lapply( 1:qa, FUN = function(i) dcorFun(x = Xb[, i], k = j, b)^2 ) )
    Vstar <- sort(Vstar)
    la <- ceiling( (1 - alpha) * qa )
    sqrt( (b - j) * Vstar[la] / (n - j) )
  }

  optimal.block <- function(x, lags, MaxLag) {
    r <- 2
    rr <- seq(-r, r, by = 1)
    bsmall <- MaxLag + 4
    bbig <- bsmall + 20
    bseq <- seq(bsmall, bbig, by = 4)
    se <- numeric( length(bseq) )
    for ( l in 1:length(bseq) ) {
      bk <- bseq[l] + rr
      Tneighb <- sapply( 1:length(bk), function(m) test(x, j = lags, b = bk[m]) )
      se[l] <- sqrt( sum( abs( Tneighb - mean(Tneighb) )^2 ) / (length(bk) - 1) )
    }
    bseq[ order(se)[1] ]
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
    fe_call <- as.call( c( list( as.name("foreach"), i = i, .combine = "c", .packages = c("dcov") ) ) )
    fe <- eval(fe_call)
    cv <- fe %dopar% test( x, i, b = optimal.block(x, lags = i, MaxLag) )
    stopCluster(cl)

  } else {
    cv <- sapply( 1:MaxLag, function(i) test( x, i, b = optimal.block(x, lags = i, MaxLag) ) )
  }

  cv
}


















# SubsCV <- function(x,MaxLag,parallel=FALSE){
#  n <- length(x)
#  if (((n-MaxLag)<0) || ((n-MaxLag)<4) || ((n-MaxLag)<=25)) stop("Give bigger sample size n")
#  dcorFun <- function(x,k,b){
#      xA <- x[1:(b-k)]
#      xB <- x[(1+k):b]
#      return(dcor(xA,xB))
#  }
#  test <- function(x,j,b){
#   q <- n-b+1
#   q <- ifelse(q<0,abs(q),q)
#   Xb <- sapply(1:q, function(i) x[i:(i+b-1)])
#   Vx <- unlist(lapply(1:q, FUN=function(i) dcorFun(x=Xb[,i],k=j,b)^2))
#   Vstar <- sort(Vx)
#   la <- ceiling((1-0.05)*q)
#   critical.value <- sqrt((b-j)*Vstar[la]/(n-j))
#   return(critical.value)
#  }
#  optimal.block <- function(x,lags,MaxLag){
#   r <- 2
#   bsmall <- MaxLag+4
#   bbig <- bsmall+20
#   bseq <- seq(bsmall,bbig,by=4)
#   se <- rep(0,length(bseq))
#   for (l in 1:length(bseq)){
#    bk <- bseq[l]+seq(-r,r,by=1)
#    Tneighb <- sapply(1:length(bk), function(m) test(x,j=lags,b=bk[m]))
#    se[l] <- (sum(abs(Tneighb-mean(Tneighb))^2)/(length(bk)-1))^(1/2)
#   }
#   smallest <- order(se)[1]
#   b.star <- bseq[smallest]
#   return(b.star)
#  }
#  if(parallel==TRUE){
#   closeAllConnections()
#   #cl <- makeCluster(detectCores())
#   cl <- makeCluster(2)
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- 1:MaxLag
#   fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="c",.export=c("dcor")) ))
#   fe <- eval(fe_call)
#   cv <- fe %dopar% test(x,i,b=optimal.block(x,lags=i,MaxLag))
#   stopCluster(cl)
#  } else {
#  cv=sapply(1:MaxLag, function(i) test(x,i,b=optimal.block(x,lags=i,MaxLag)))
#  }
#  return(cv)
# }
