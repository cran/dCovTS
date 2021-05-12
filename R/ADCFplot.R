ADCFplot <- function(x, MaxLag = 15, alpha = 0.05, b = 499, bootMethod =
            c("Wild Bootstrap", "Subsampling", "Independent Bootstrap"), ylim = NULL, main = NULL) {

  if ( b <= 0 )  stop( "No plot is given for b<=0." )
  if ( MaxLag == 0 )  stop( "MaxLag must be greater than 0." )
  series <- deparse( substitute(x) )
  bootMethod <- match.arg(bootMethod)
  if ( missing(bootMethod) )  bootMethod <- "Wild Bootstrap"

  n <- length(x)
  adcor <- ADCF(x, MaxLag, unbiased = FALSE)

  if ( bootMethod == "Wild Bootstrap" ) {
    cv <- RbootCV(n, MaxLag, alpha, b, parallel = TRUE)
  } else if ( bootMethod == "Independent Bootstrap" ) {
    cv <- OrdinaryBootCV(n, MaxLag, alpha, b, parallel = TRUE)
  } else {
    if ( ( (n - MaxLag) < 0 )  ||  ( (n - MaxLag) < 4 ) || ( (n - MaxLag) <= 25 ) )  stop( "Give bigger sample size n." )
    cv <- SubsCV(x, MaxLag, alpha, parallel = TRUE)
  }

  r1 <- max(cv, 1)
  if ( is.null(ylim) ) ylim <- c(0, r1)

  if ( length(cv) == 1 ) {
    if ( is.null(main) ) {
      plot(0:MaxLag, adcor, type = "n", main = paste("Series", series), xlab = "Lag",
           ylab = "ADCF", ylim = ylim, cex.lab = 1.2, cex.axis = 1.2)
    } else  plot(0:MaxLag, adcor, type = "n", main = main, xlab = "Lag", ylab = "ADCF",
                  ylim = ylim, cex.lab = 1.2, cex.axis = 1.2)
    for ( i in seq(0, MaxLag, by = 1) )  segments(i, 0, i, adcor[i + 1])
    points(0:MaxLag, rep(cv, MaxLag + 1), type = "l", lty = 3, lwd = 2, col = "blue")
  } else {
    if ( is.null(main) ) {
      plot(1:MaxLag, adcor[-1], type = "n", main = paste("Series", series), xlab = "Lag",
            ylab = "ADCF", ylim = ylim, cex.lab = 1.2, cex.axis = 1.2)
    } else  plot(1:MaxLag, adcor[-1], type = "n", main = main, xlab = "Lag", ylab = "ADCF",
                 ylim = ylim, cex.lab = 1.2, cex.axis = 1.2)
    for ( i in seq(1, MaxLag, by = 1) )  segments(i, 0, i, adcor[i + 1])
    points(1:MaxLag, cv, type = "l", lty = 3, lwd = 2, col = "blue")
  }

  result <- list(ADCF = adcor, bootMethod = bootMethod, critical.values = cv)
  return(result)
}


















# ADCFplot <- function (x, MaxLag = 15, ylim=NULL, main = NULL, bootMethod = c("Wild Bootstrap",
#     "Subsampling","Independent Bootstrap"), b = 499)
# {
#     if (b <= 0)
#         stop("No plot is given for b<=0")
#     if (MaxLag==0)
#         stop("MaxLag must be greater than 0")
#     series <- deparse(substitute(x))
#     bootMethod <- match.arg(bootMethod)
#     if (missing(bootMethod))
#         bootMethod = "Wild Bootstrap"
#     n <- length(x)
#     adcor <- ADCF(x, MaxLag, unbiased=FALSE)
#     if (bootMethod == "Wild Bootstrap") {
#         cv <- RbootCV(n, MaxLag, b = b, parallel = TRUE)
#     }
#     else if (bootMethod == "Independent Bootstrap") {
#         cv <- OrdinaryBootCV(n, MaxLag, b = b, parallel = TRUE)
#     }
#     else {
#         if (((n - MaxLag) < 0) || ((n - MaxLag) < 4) || ((n -
#             MaxLag) <= 25))
#             stop("Give bigger sample size n")
#         cv <- SubsCV(x,MaxLag,parallel=TRUE)
#     }
#     r1 <- max(cv, 1)
#     if (is.null(ylim)) ylim=c(0,r1)
#     if (length(cv) == 1) {
#         if (is.null(main)) {
#             plot(0:MaxLag, adcor, type = "n", main = paste("Series",
#                 series), xlab = "Lag", ylab = "ADCF", ylim = ylim)
#         }
#         else {
#             plot(0:MaxLag, adcor, type = "n", main = main, xlab = "Lag",
#                 ylab = "ADCF", ylim = ylim)
#         }
#         for (i in seq(0, MaxLag, by = 1)) {
#             segments(i, 0, i, adcor[i + 1])
#         }
#         points(0:MaxLag, rep(cv, MaxLag + 1), type = "l", lty = 3,
#             lwd = 2, col = "blue")
#     }
#     else {
#         if (is.null(main)) {
#             plot(1:MaxLag, adcor[-1], type = "n", main = paste("Series",
#                 series), xlab = "Lag", ylab = "ADCF", ylim = ylim)
#         }
#         else {
#             plot(1:MaxLag, adcor[-1], type = "n", main = main,
#                 xlab = "Lag", ylab = "ADCF", ylim = ylim)
#         }
#         for (i in seq(1, MaxLag, by = 1)) {
#             segments(i, 0, i, adcor[i + 1])
#         }
#         points(1:MaxLag, cv, type = "l", lty = 3, lwd = 2, col = "blue")
#     }
#     result <- list(ADCF = adcor, bootMethod = bootMethod, critical.values = cv)
#     return(result)
# }
