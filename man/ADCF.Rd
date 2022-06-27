\name{ADCF}
\alias{ADCF}
\title{
Auto-Distance Correlation Function
}

\description{
Computes the auto-distance correlation function of a univariate time series. It also computes the bias-corrected
estimator of (squared) auto-distance correlation.
}
\usage{
ADCF(x, MaxLag = 15, unbiased = FALSE)
}

\arguments{
\item{x}{
A numeric vector or univariate time series.
}
\item{MaxLag}{
Maximum lag order at which to calculate the \code{ADCF}. Default is 15.
}
\item{unbiased}{
Logical value. If unbiased = TRUE, the bias-corrected estimator of squared auto-distance correlation is returned.
Default value is FALSE.
}
}
\details{
Distance covariance and correlation firstly introduced by Szekely et al. (2007) are new measures of dependence between two
random vectors. Zhou (2012) extended this measure to univariate time series framework.

For a univariate time series, \code{\link{ADCF}} computes the auto-distance correlation function, \eqn{R_X(j)},
between \eqn{\{X_t\}} and \eqn{\{X_{t+j}\}}, whereas \code{\link{ADCV}} computes the auto-distance covariance function
between them, denoted by \eqn{V_X(j)}. Formal definition of \eqn{R_X(\cdot)} and \eqn{V_X(\cdot)} can be found in Zhou (2012)
and Fokianos and Pitsillou (2017). The empirical auto-distance correlation function, \eqn{\hat{R}_X(j)}, is computed as the
positive square root of \deqn{ \hat{R}_X^2(j)=\frac{\hat{V}_X^2(j)}{\hat{V}_X^2(0)}, \quad j=0, \pm 1, \pm 2, \dots}%{ \hat{R}_X^2(j) = \hat{V}_X^2(j) / \hat{V}_X^2(0),    j=0, \pm 1, \pm 2, \dots}

for \eqn{\hat{V}_X^2(0) \neq 0} and zero otherwise, where \eqn{\hat{V}_X(\cdot)} is a function of the double
centered Euclidean distance matrices of the sample \eqn{X_t} and its lagged sample \eqn{X_{t+j}} (see \code{\link{ADCV}}
for more details). Theoretical properties of this measure can be found in Fokianos and Pitsillou (2017).

If unbiased = TRUE, \code{\link{ADCF}} computes the bias-corrected estimator of the squared auto-distance correlation,
\eqn{\tilde{R}_X^2(j)}, based on the unbiased estimator of auto-distance covariance function, \eqn{\tilde{V}_X^2(j)}.
The definition of \eqn{\tilde{V}_X^2(j)} relies on the U-centered matrices proposed by Szekely and Rizzo (2014)
(see \code{\link{ADCV}} for a brief description).

\code{\link{mADCF}} computes the auto-distance correlation function of a multivariate time series.
}
\value{
Returns a vector, whose length is determined by \code{MaxLag}, and contains the biased estimator of \code{ADCF}
or the bias-corrected estimator of squared \code{ADCF}.
}
\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature Review of Distance Correlation and Its
Applications to Time Series. \emph{International Statistical Review}, 87, 237-262.

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise dependence in time series.
\emph{Technometrics}, 159(2), 262-3270.

Huo, X. and G. J. Szekely. (2016). Fast Computing for Distance Covariance.
\emph{Technometrics}, 58, 435-447.

Pitsillou M. and Fokianos K. (2016). dCovTS: Distance Covariance/Correlation for Time Series.
\emph{R Journal}, 8, 324-340.

Szekely, G. J. and M. L. Rizzo (2014). Partial distance correlation with methods for dissimilarities.
\emph{The Annals of Statistics}, 42, 2382-2412.

Szekely, G. J. and M. L. Rizzo and N. K. Bakirov (2007). Measuring and testing dependence by correlation
of distances. \emph{The Annals of Statistics}, 35, 2769-2794.

Zhou, Z. (2012). Measuring nonlinear dependence in time series, a distance correlation approach.
\emph{Journal of Time Series Analysis} 33, 438-457.
}

\author{
Maria Pitsillou, Michail Tsagris and Konstantinos Fokianos.
}
\note{
Based on the definition of \code{ADCF}, one can observe that \eqn{R^2_X(j)=R^2_X(-j)} \eqn{\forall} \eqn{j},
and so results based on negative lags are omitted.
}

\seealso{
\code{\link{ADCFplot}}, \code{\link{ADCV}}, \code{\link{mADCF}}
}
\examples{
x <- rnorm(100)
ADCF(x)

ADCF(ldeaths, 18)

ADCF(mdeaths, unbiased = TRUE)
}
\keyword{ ts }
