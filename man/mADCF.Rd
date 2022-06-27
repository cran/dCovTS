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
Multivariate time series.
}
\item{lags}{
The lag order at which to calculate the \code{mADCF}. No default is given.
This can be a single number or a vector of numbers with different lag orders.
}
\item{unbiased}{
A logical value. If unbiased = TRUE, the individual elements of auto-distance
correlation matrix correspond to the bias-corrected estimators of squared
auto-distance correlation functions. Default value is FALSE.
}
\item{output}{
A logical value. If output=FALSE, no output is given. Default value is TRUE.
}
}
\details{
If \eqn{\textbf{X}_t=(X_{t;1}, \dots, X_{t;d})'} is a multivariate time series
of dimension \eqn{d}, then \code{\link{mADCF}} computes the sample auto-distance
correlation matrix, \eqn{\hat{R}(\cdot)}, of \eqn{\textbf{X}_t}. It is defined by
\deqn{ \hat{R}(j) = [\hat{R}_{rm}(j)]_{r,m=1}^d, \quad j=0, \pm 1, \pm 2, \dots,
}

where \eqn{\hat{R}_{rm}(j)} is the biased estimator of the so-called pairwise
auto-distance correlation function between
\eqn{X_{t;r}} and \eqn{X_{t+j;m}} given by the positive square root of
\deqn{ \hat{R}_{rm}^2(j) = \frac{\hat{V}_{rm}^2(j)}{\hat{V}_{rr}(0)\hat{V}_{mm}(0)}
}{%
\hat{R}_{rm}^2(j) = \hat{V}_{rm}^2(j) / (\hat{V}_{rr}(0)\hat{V}_{mm}(0))
}
for \eqn{\hat{V}_{rr}(0)\hat{V}_{mm}(0) \neq 0} and zero otherwise.

\eqn{\hat{V}_{rm}(j)} is the \eqn{(r,m)} element of the corresponding
\code{\link{mADCV}} matrix at lag \eqn{j}.
Formal definition and more details can be found in Fokianos and Pitsillou (2017).

If unbiased = TRUE, \code{\link{mADCF}} returns a matrix that contains the
bias-corrected estimators of squared pairwise
auto-distance correlation functions.
}
\value{
If \code{lags} is a single number then the function will return a matrix.
If \code{lags} is a vector of many values the function will return an array.
For either case, the matrix (matrices) will contain either the biased
estimators of the pairwise auto-distance correlation functions or the
bias-corrected estimators of squared pairwise auto-distance correlation
functions at lag, \eqn{j}, determined by the argument \code{lags}.
}
\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature
Review of Distance Correlation and Its Applications to Time Series.
\emph{International Statistical Review}, 87, 237-262.

Fokianos K. and Pitsillou M. (2018). Testing independence for multivariate
time series via the auto-distance correlation matrix.
\emph{Biometrika}, 105, 337-352.

Huo, X. and G. J. Szekely. (2016). Fast Computing for Distance Covariance.
\emph{Technometrics}, 58, 435-447.

Pitsillou M. and Fokianos K. (2016). dCovTS: Distance Covariance/Correlation
for Time Series. \emph{R Journal}, 8, 324-340.
}
\author{
Maria Pitsillou, Michail Tsagris and Konstantinos Fokianos.
}
\seealso{
\code{ \link{ADCF}}, \code{\link{mADCV} }
}
\examples{
x <- matrix( rnorm(200), ncol = 2 )

mADCF(x, 2)

mADCF(x, -2)

mADCF(x, lags = 4, unbiased = TRUE)
}
\keyword{ ts }
