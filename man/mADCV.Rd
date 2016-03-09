\name{mADCV}
\alias{mADCV}
\title{
Distance Covariance Matrix
}
\description{
Computes the sample distance covariance matrices of a multivariate time series.
}
\usage{
mADCV(x, lags, output=TRUE)
}
\arguments{
  \item{x}{
multivariate time series.
}
  \item{lags}{
lag order at which to calculate the \code{mADCV}. No default is given.
}
  \item{output}{
logical value. If output=FALSE, no output is given. Default value is TRUE.
}
}
\details{
Suppose that \eqn{\textbf{X}_t=(X_{t;1}, \dots, X_{t;d})'} is a multivariate time series of dimension \eqn{d}. Then, 
\code{mADCV} computes the \eqn{d \times d} sample distance covariance matrix, \eqn{\hat{V}(\cdot)}, of \eqn{\textbf{X}_t} given by
\deqn{ \hat{V}(j) = [\hat{V}_{rm}(j)]_{r,m=1}^d , j~~=~~0, \pm 1, \pm 2, \dots
}
where \eqn{\hat{V}_{rm}(j)} denotes the pairwise sample auto-distance covariance function between \eqn{X_{t;r}} and \eqn{X_{t-|j|;m}}
whose definition is given analogously as in the univariate case \code{\link{ADCV}}. Formal definitions and theoretical properties of \code{mADCV}
can be found in Fokianos and Pitsillou (2016).
}
\value{
Returns the sample auto-distance covariance matrix at lag, \eqn{j}, determined by the argument \code{lags}.
}
\references{
Fokianos K. and M. Pitsillou (2016). On multivariate auto-distance covariance and correlation functions. Submitted for publication.
}
\author{
Maria Pitsillou and Konstantinos Fokianos
}


\seealso{
\code{\link{ADCV}},\code{\link{mADCF}}
}
\examples{
x <- MASS::mvrnorm(100,rep(0,2),diag(2))
mADCV(x,lags=1)
mADCV(x,lags=15)

y <- as.ts(swiss)
mADCV(y,15)
}


\keyword{ ts }