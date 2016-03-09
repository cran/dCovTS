\name{ibmSp500}
\alias{ibmSp500}
\docType{data}
\title{
Monthly returns of IBM and S&P 500 composite index
}
\description{
The monthly returns of the stocks of International Business Machines (IBM) and the S&P 500 composite index from January
1926 to December 2011.
}
\usage{ibmSp500}
\format{
  A data frame with 1032 observations on the following 3 variables.
  \describe{
    \item{\code{date}}{a numeric vector}
    \item{\code{ibm}}{a numeric vector}
    \item{\code{sp}}{a numeric vector}
  }
}
\source{
The data is a combination of two datasets:
\itemize{
\item The first 612 observations are in Tsay (2010) (see \url{http://faculty.chicagobooth.edu/ruey.tsay/teaching/fts3/})
\item The rest 420 observations are in Tsay (2014) (see \url{http://faculty.chicagobooth.edu/ruey.tsay/teaching/mtsbk/})
}
}                                                                                                                                            
\references{
Tsay, R. S. (2010). \emph{Analysis of Financial Time Series}. Hoboken, NJ: Wiley. Third edition.

Tsay, R. S. (2014). \emph{Multivariate Time Series Analysis with R and Financial Applications}. Hoboken, NJ:
Wiley.
}
\examples{
attach(ibmSp500)

series <- tail(ibmSp500[,2:3],400)
lseries <- log(series+1)
\dontrun{
mADCFplot(lseries,MaxLag=18)
mADCFplot(lseries^2,MaxLag=18)
acf(lseries,lag.max=18)
acf(lseries^2,lag.max=18)}

}
\keyword{datasets}
