mRbootCV <- function(n, qa, MaxLag, alpha, b, parallel = FALSE) {

  x <- Rfast::matrnorm(n, qa)
  if ( missing(MaxLag)  ||  MaxLag < 0 )  stop( "'MaxLag' must be greater than 0." )
  A0 <- crossDist(x, 0)$A
  Atilde0 <- lapply(1:qa, FUN = function(j)  ATilde( A0[[ j ]] ) )

  rstar <- function(k) {
    cross <- crossDist(x, k)
    A <- cross$A
    B <- cross$B
    Atilde <- Btilde <- list()
    for (i in 1:qa)  {
      Atilde[[ i ]] <- ATilde( A[[ i ]] )
    	Btilde[[ i ]] <- ATilde( B[[ i ]] )
    }

    ela <- expand.grid(1:qa, 1:qa)
    dvarx <- numeric( dim(ela)[1] )
    for ( i in 1:dim(ela)[1] ) {
      dvarx[i] <- mean( Atilde0[[ ela[i, 1] ]]^2 ) * mean( Atilde0[[ ela[i, 2] ]]^2 )
    }
    dvarx <- dvarx^0.25

    Rarray <- function(Atilde, Btilde, Atilde0, k, dvarx) {
      ela <- expand.grid(1:qa, 1:qa)
      Rm <- matrix(nrow = qa, ncol = qa)
      Wtstar <- Rfast::Rnorm(n - k)
      down <- (n - k) * dvarx
      for ( i in 1:dim(ela)[1] ) {
        com <- Atilde[[ ela[i, 1] ]] * Btilde[[ ela[i, 2] ]]
        dcova <- sum( Wtstar %*% com * Wtstar )
        Rm[ela[i, 1], ela[i, 2]] <- sqrt( dcova ) / down[i]
      }
      Rm
    }

    result <- replicate( b, Rarray(Atilde, Btilde, Atilde0, k, dvarx) )

    cv <- numeric( dim(ela)[1] )
    for ( i in 1:dim(ela)[1] ) {
      quant <- quantile( result[ela[i, 1], ela[i, 2], ], 1 - alpha )
      pv <- mean( result[ela[i, 1], ela[i, 2], ] >= quant )
      pvadj <- p.adjust( pv, method = "fdr" )
      cv[i] <- quantile( result[ela[i, 1], ela[i, 2], ], 1 - pvadj )
    }
    cv
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
    fe_call <- as.call( c( list( as.name("foreach"), i = i, .combine = "rbind",
                        .export = c("crossDist", "ATilde"), .packages = "Rfast" ) ) )
    fe <- eval(fe_call)
    Rstar <- fe %dopar% rstar(i)
    stopCluster(cl)
    mcv <- matrix( Rfast::colMaxs(Rstar, TRUE), nrow = qa, ncol = qa, byrow = TRUE)

  } else {
    res <- lapply( 1:MaxLag,FUN = function(i) rstar(i) )
    Rstar <- t( sapply(1:2, FUN = function(i) unlist( res[[ i ]] ) ) )
    mcv <- matrix( sapply(1:qa^2, function(i) max( Rstar[[ i ]]) ), nrow = qa, ncol = qa, byrow = TRUE )
  }
  mcv
}

















# mRbootCV <- function(n,q,MaxLag,b,parallel=FALSE){
#  x <- replicate(q,rnorm(n))
#  if (missing(MaxLag) || MaxLag < 0)
#       stop("'MaxLag' must be greater than 0")
#  A0 <- crossDist(x,0)$A
#  Atilde0 <- lapply(1:q, FUN=function(j) ATilde(A0[[j]]))
#   rstar <- function(k){
#    A <- crossDist(x,k)$A
#    B <- crossDist(x,k)$B
#    Atilde <- lapply(1:q, FUN=function(j) ATilde(A[[j]]))
#    Btilde <- lapply(1:q, FUN=function(j) ATilde(B[[j]]))
#    cv <- vector()
#    Rarray <- function(Atilde,Btilde,Atilde0,k){
#     Wtstar <- rbind(rnorm(n-k))
#     Wt <- rbind(rnorm(n))
#     Rm <- matrix(NA,q,q)
#      for (i in 1:q){
#      for (j in 1:q){
#        dcov <- sqrt((Wtstar%*%(Atilde[[i]]*Btilde[[j]])%*%t(Wtstar))/((n-k)^2))
#        dvarx <- sqrt(mean((Atilde0[[i]]*Atilde0[[i]]))*mean((Atilde0[[j]]*Atilde0[[j]])))
#        Rm[i,j] <- dcov/sqrt(dvarx)  
#       }
#      }
#     return(Rm)
#    }
#    result <- replicate(b,Rarray(Atilde,Btilde,Atilde0,k))
#    s <- 1
#    for(i in 1:q){
#     for (j in 1:q){
#      quant <- quantile(result[i,j,],0.95)
#      pv <- mean(result[i,j,]>=quant)
#      pvadj <- p.adjust(pv,method="fdr")
#      cv[s] <- quantile(result[i,j,],1-pvadj)
#      s <- s+1
#     }
#    }
#    return(cv)
#   }
#  if(parallel==TRUE){
#   closeAllConnections()
#   #cl <- makeCluster(detectCores())
#   cl <- makeCluster(2)
#   registerDoParallel(cl)
#   clusterSetRNGStream(cl = cl, iseed = 9182)
#   i <- 1:MaxLag
#   fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="rbind",.export=c("crossDist","ATilde")) ))
#   fe <- eval(fe_call)
#   Rstar <- fe %dopar% rstar(i)
#   stopCluster(cl)
#   mcv <- matrix(sapply(1:q^2, function(i) max(Rstar[,i])),nrow=q,ncol=q,byrow=T)
#  }
#  else {
#  res <- lapply(1:MaxLag,FUN=function(i) rstar(i))
#  Rstar <- t(sapply(1:2,FUN=function(i) unlist(res[[i]])))
#  mcv <- matrix(sapply(1:q^2, function(i) max(Rstar[[i]])),nrow=q,ncol=q,byrow=T)
#  }
#  return(mcv)
# }
