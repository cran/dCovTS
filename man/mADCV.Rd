\name{mADCV}
\alias{mADCV}
\title{
Auto-Distance Covariance Matrix
}
\description{
Computes the sample auto-distance covariance matrices of a multivariate time series.
}
\usage{
mADCV(x, lags, unbiased = FALSE, output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
multivariate time series.
}
  \item{lags}{
lag order at which to calculate the \code{mADCV}. No default is given.
}
  \item{unbiased}{
logical value. If unbiased = TRUE, the individual elements of auto-distance covariance matrix correspond to
 the unbiased estimators of squared auto-distance covariance functions. Default value is FALSE.
}
  \item{output}{
logical value. If output=FALSE, no output is given. Default value is TRUE.
}
}
\details{
Suppose that \eqn{\textbf{X}_t=(X_{t;1}, \dots, X_{t;d})'} is a multivariate time series of dimension \eqn{d}. Then, 
\code{\link{mADCV}} computes the \eqn{d \times d} sample auto-distance covariance matrix, \eqn{\hat{V}(\cdot)}, of \eqn{\textbf{X}_t} given by
\deqn{ \hat{V}(j) = [\hat{V}_{rm}(j)]_{r,m=1}^d , \quad j=0, \pm 1, \pm 2, \dots,
}
where \eqn{\hat{V}_{rm}(j)} denotes the biased estimator of the pairwise auto-distance covariance function between \eqn{X_{t;r}} and \eqn{X_{t+j;m}}.
The definition of \eqn{\hat{V}_{rm}(j)} is given analogously as in the univariate case (see \code{\link{ADCV}}). 
Formal definitions and theoretical properties of auto-distance covariance matrix can be found in Fokianos and Pitsillou (2016).

If unbiased = TRUE, \code{\link{mADCV}} computes the matrix, \eqn{\tilde{V}^{(2)}(j)}, whose elements correspond to the unbiased estimators of
squared pairwise auto-distance covariance functions, namely
\deqn{ \tilde{V}^{(2)}(j) = [\tilde{V}^2_{rm}(j)]_{r,m=1}^d , \quad j=0, \pm 1, \pm 2, \dots.
}
The definition of \eqn{\tilde{V}_{rm}^2(\cdot)} is defined analogously as explained in the univariate case (see \code{\link{ADCV}}).
}
\value{
Returns a matrix containing either the biased estimators of the pairwise auto-distance covariance functions
or the unbiased estimators of squared pairwise auto-distance covariance functions
at lag, \eqn{j}, determined by the argument \code{lags}.
}
\references{
Fokianos K. and M. Pitsillou (2016). Testing pairwise independence for multivariate
time series by the auto-distance correlation matrix. Submitted for publication.
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
mADCV(y,15,unbiased=TRUE)
}


\keyword{ ts }