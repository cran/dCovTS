\name{mADCF}
\alias{mADCF}
\title{
Distance Correlation Matrix 
}
\description{
Computes the distance correlation matrix of a multivariate time series.
}
\usage{
mADCF(x, lags, output=TRUE)
}

\arguments{
  \item{x}{
multivariate time series.
}
  \item{lags}{
lag order at which to calculate the \code{mADCF}. No default is given.
}
  \item{output}{
logical value. If output=FALSE, no output is given. Default value is TRUE.
}
}
\details{
If \eqn{\textbf{X}_t=(X_{t;1}, \dots, X_{t;d})'} is a multivariate time series of dimension \eqn{d}, then \code{mADCF}
computes the sample distance correlation matrix, \eqn{\hat{R}(\cdot)}, of \eqn{\textbf{X}_t}. It is defined by
\deqn{ \hat{R}(j) = [\hat{R}_{rm}(j)]_{r,m=1}^d  ~~~~, j~~=~~0, \pm 1, \pm 2, \dots  
}
where \eqn{\hat{R}_{rm}(j)} is the sample auto-distance correlation function between \eqn{X_{t;r}} and \eqn{X_{t-|j|;m}} given as
the positive square root of 
\deqn{ \hat{R}_{rm}^2(j) = \frac{\hat{V}_{rm}^2(j)}{\hat{V}_{rr}(0)\hat{V}_{mm}(0)}  
}
for \eqn{\hat{V}_{rr}(0)\hat{V}_{mm}(0) \neq 0} and zero otherwise.

\eqn{\hat{V}_{rm}(j)} is the \eqn{(r,m)}-element of the corresponding \code{\link{mADCV}} matrix at lag \eqn{j}. Formal definition and 
more details can be found in Fokianos and Pitsillou (2016).
}
\value{
Returns the sample auto-distance correlation matrix at lag, \eqn{j}, determined by the argument \code{lags}.
}
\references{
Fokianos K. and M. Pitsillou (2016). On multivariate auto-distance covariance and correlation functions. Submitted for publication.
}
\author{
Maria Pitsillou and Konstantinos Fokianos
}

\seealso{
\code{\link{ADCF}},\code{\link{mADCV}}
}
\examples{
x <- MASS::mvrnorm(100,rep(0,2),diag(2))
mADCF(x,2)

mADCF(x,-2)
}
\keyword{ ts }
