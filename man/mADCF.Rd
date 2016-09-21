\name{mADCF}
\alias{mADCF}
\title{
Auto-Distance Correlation Matrix 
}
\description{
Computes the auto-distance correlation matrix of a multivariate time series.
}
\usage{
mADCF(x, lags, unbiased = FALSE, output = TRUE)
}

\arguments{
  \item{x}{
multivariate time series.
}
  \item{lags}{
lag order at which to calculate the \code{mADCF}. No default is given.
}
  \item{unbiased}{
logical value. If unbiased = TRUE, the individual elements of auto-distance correlation matrix correspond to
 the bias-corrected estimators of squared auto-distance correlation functions. Default value is FALSE.
}
  \item{output}{
logical value. If output=FALSE, no output is given. Default value is TRUE.
}
}
\details{
If \eqn{\textbf{X}_t=(X_{t;1}, \dots, X_{t;d})'} is a multivariate time series of dimension \eqn{d}, then \code{\link{mADCF}}
computes the sample auto-distance correlation matrix, \eqn{\hat{R}(\cdot)}, of \eqn{\textbf{X}_t}. It is defined by
\deqn{ \hat{R}(j) = [\hat{R}_{rm}(j)]_{r,m=1}^d, \quad j=0, \pm 1, \pm 2, \dots,  
}

where \eqn{\hat{R}_{rm}(j)} is the biased estimator of the so-called pairwise auto-distance correlation function between 
\eqn{X_{t;r}} and \eqn{X_{t+j;m}} given by the positive square root of 
\deqn{ \hat{R}_{rm}^2(j) = \frac{\hat{V}_{rm}^2(j)}{\hat{V}_{rr}(0)\hat{V}_{mm}(0)}  
}{%
\hat{R}_{rm}^2(j) = \hat{V}_{rm}^2(j) / (\hat{V}_{rr}(0)\hat{V}_{mm}(0))
}
for \eqn{\hat{V}_{rr}(0)\hat{V}_{mm}(0) \neq 0} and zero otherwise.

\eqn{\hat{V}_{rm}(j)} is the \eqn{(r,m)} element of the corresponding \code{\link{mADCV}} matrix at lag \eqn{j}. Formal definition and 
more details can be found in Fokianos and Pitsillou (2016).

If unbiased = TRUE, \code{\link{mADCF}} returns a matrix that contains the bias-corrected estimators of squared pairwise auto-distance correlation
functions, namely
\deqn{ \tilde{R}^{(2)}(j) = [\tilde{R}^2_{rm}(j)]_{r,m=1}^d, \quad j=0, \pm 1, \pm 2, \dots.
}
\eqn{\tilde{R}_{rm}^2(\cdot)} are based on the unbiased estimator of pairwise auto-distance covariance, \eqn{\tilde{V}_{rm}^2(\cdot)}. 
The definition of \eqn{\tilde{V}_{rm}^2(\cdot)} can be found in \code{\link{mADCV}}.
}
\value{
Returns a matrix containing either the biased estimators of the pairwise auto-distance correlation functions
or the bias-corrected estimators of squared pairwise auto-distance correlation functions
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
\code{\link{ADCF}},\code{\link{mADCV}}
}
\examples{
x <- MASS::mvrnorm(100,rep(0,2),diag(2))
mADCF(x,2)

mADCF(x,-2)

mADCF(x,lags=4,unbiased=TRUE)
}
\keyword{ ts }
