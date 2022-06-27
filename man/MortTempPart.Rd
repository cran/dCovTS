\name{MortTempPart}
\alias{MortTempPart}
\docType{data}
\title{
Cardiovascular mortality, temperature and pollution data in Los Angeles County
}
\description{Cardiovascular mortality data measured daily in Los Angeles
County over the 10 year period 1970-1979. Temperature series and pollutant
particulate series corresponding to mortality data are also given.
}
\usage{MortTempPart}
\format{
A data frame with 508 observations on the following 3 variables.
\describe{
\item{\code{cmort}}{A numeric vector.}
\item{\code{tempr}}{A numeric vector.}
\item{\code{part}}{A numeric vector.}
}
}
\references{
Shumway R. H. and D. S. Stoffer (2011).
\emph{Time Series Analysis and Its Applications With R Examples}.
New York: Springer. Third Edition.
}
%\note{
%Due to restrictions in the use of parallel on CRAN, we removed the example.
%An example function for this dataset is "mADCFplot(MortTempPart[1:100, ], MaxLag = 10)".
%}
\examples{
\donttest{
### mADCFplot(MortTempPart[1:100, ], MaxLag = 10)
}
}
\keyword{datasets}
