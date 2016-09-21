\name{ADCF}
\alias{ADCF}

\title{
Auto-Distance Correlation Function
}
\description{
Computes the auto-distance correlation function of a univariate time series. It also computes the bias-corrected estimator of (squared) auto-distance correlation.
}
\usage{
ADCF(x, MaxLag = 15, unbiased = FALSE)
}

\arguments{
  \item{x}{
numeric vector or univariate time series.
}
  \item{MaxLag}{
maximum lag order at which to calculate the \code{ADCF}. Default is 15. 
}
  \item{unbiased}{
logical value. If unbiased = TRUE, the bias-corrected estimator of squared auto-distance correlation is returned. Default value is FALSE.
}
}
\details{
Distance covariance and correlation firstly introduced by Szekely et al. (2007) are new measures of dependence between two random vectors. Zhou (2012) extended this measure
to time series framework.

For a univariate time series, \code{\link{ADCF}} computes the auto-distance correlation function, \eqn{R_X(j)}, between \eqn{\{X_t\}} and \eqn{\{X_{t+j}\}}, 
whereas \code{\link{ADCV}} computes the auto-distance covariance function between them, denoted by \eqn{V_X(j)}. Formal definition
of \eqn{R_X(\cdot)} and \eqn{V_X(\cdot)} can be found in Zhou (2012) and Fokianos and Pitsillou (2016). The empirical auto-distance correlation function, \eqn{\hat{R}_X(j)}, is computed as the positive
square root of
\deqn{ \hat{R}_X^2(j)=\frac{\hat{V}_X^2(j)}{\hat{V}_X^2(0)}, \quad j=0, \pm 1, \pm 2, \dots}{%
\hat{R}_X^2(j) = \hat{V}_X^2(j) / \hat{V}_X^2(0),    j=0, \pm 1, \pm 2, \dots} 

 for \eqn{\hat{V}_X^2(0) \neq 0} and zero otherwise, where \eqn{\hat{V}_X(\cdot)} is a function of the double
centered Euclidean distance matrices of the sample \eqn{X_t} and its lagged sample \eqn{X_{t+j}} (see \code{\link{ADCV}} for more details).
Theoretical properties of this measure can be found in Fokianos and Pitsillou (2016).

If unbiased = TRUE, \code{\link{ADCF}} computes the bias-corrected estimator of the squared auto-distance correlation, \eqn{\tilde{R}_X^2(j)},
based on the unbiased estimator of auto-distance covariance function, \eqn{\tilde{V}_X^2(j)}. The definition of \eqn{\tilde{V}_X^2(j)} relies on the U-centered matrices proposed by 
Szekely and Rizzo (2014) (see \code{\link{ADCV}} for a brief description).

\code{\link{mADCF}} computes the auto-distance correlation function of a multivariate time series.
}
\value{
Returns a vector, whose length is determined by \code{MaxLag}, and contains the biased estimator of \code{ADCF} or the bias-corrected
estimator of squared \code{ADCF}. 
}
\references{
Fokianos K. and M. Pitsillou (2016). Consistent testing for pairwise dependence in time series. \emph{Technometrics},
\url{http://dx.doi.org/10.1080/00401706.2016.1156024}.

Szekely, G. J. and M. L. Rizzo (2014). Partial distance correlation with methods for dissimilarities. 
\emph{The Annals of Statistics} \eqn{\textbf{42}}{42}, 2382-2412,
\url{dx.doi.org/10.1214/14-AOS1255}.

Szekely, G. J. and M. L. Rizzo and N. K. Bakirov (2007). Measuring and testing dependence 
by correlation of distances. \emph{The Annals of Statistics} \eqn{\textbf{35}}{35}, 2769-2794,
\url{http://dx.doi.org/10.1214/009053607000000505}.

Zhou, Z. (2012). Measuring nonlinear dependence in time series, a distance correlation approach.
\emph{Journal of Time Series Analysis} \eqn{\textbf{33}}{33}, 438-457,
\url{http://dx.doi.org/10.1111/j.1467-9892.2011.00780.x}.

}
\author{
Maria Pitsillou and Konstantinos Fokianos
}
\note{
Based on the definition of \code{ADCF}, one can observe that \eqn{R^2_X(j)=R^2_X(-j)} \eqn{\forall} \eqn{j}, 
and so results based on negative lags are omitted.
}

\seealso{
\code{\link{ADCFplot}}, \code{\link{ADCV}}, \code{\link{mADCF}}
}
\examples{
x <- rnorm(1000)
\dontrun{ADCF(x)}

ADCF(ldeaths,18)

ADCF(mdeaths,unbiased=TRUE)
}                                   
\keyword{ ts }
