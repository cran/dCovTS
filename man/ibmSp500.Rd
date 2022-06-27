\name{ibmSp500}
\alias{ibmSp500}
\docType{data}
\title{
Monthly returns of IBM and S&P 500 composite index
}
\description{
The monthly returns of the stocks of International Business Machines (IBM)
and the S&P 500 composite index from January 1926 to December 2011.
}
\usage{ibmSp500}
\format{
A data frame with 1,032 observations on the following 3 variables.
\describe{
\item{\code{date}}{a numeric vector}
\item{\code{ibm}}{a numeric vector}
\item{\code{sp}}{a numeric vector}
}
}
\source{
The data is a combination of two datasets:
\itemize{
\item The first 612 observations are in Tsay (2010).
\item The rest 420 observations are in Tsay (2014).
}
}
\references{
Tsay, R. S. (2010). \emph{Analysis of Financial Time Series}.
Hoboken, NJ: Wiley. Third edition.

Tsay, R. S. (2014). \emph{Multivariate Time Series Analysis with R and
Financial Applications}. Hoboken, NJ: Wiley.
}
%\note{
%Due to restrictions in the use of parallel on CRAN, we removed the example.
%Two example functions for this dataset are
%"attach(ibmSp500)" ;
%"series <- tail(ibmSp500[, 2:3], 400)" ;
%"lseries <- log(series + 1)" ;
%"mADCFplot(lseries, MaxLag = 12)" ;
%"mADCFplot(lseries^2, MaxLag = 12)".
%}
\examples{
\donttest{
### attach(ibmSp500)
### series <- tail(ibmSp500[, 2:3], 400)
### lseries <- log(series + 1)
### mADCFplot(lseries, MaxLag = 12)
### mADCFplot(lseries^2, MaxLag = 12)
}
}
\keyword{datasets}
