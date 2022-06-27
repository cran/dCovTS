\name{ADCV}
\alias{ADCV}

\title{
Auto-distance Covariance Function
}
\description{
Computes the auto-distance covariance function of a univariate time series.
It also computes the unbiased estimator of squared auto-distance covariance.
}
\usage{
ADCV(x, MaxLag = 15, unbiased = FALSE)
}
\arguments{
\item{x}{
A numeric vector or univariate time series.
}
\item{MaxLag}{
The maximum lag order at which to calculate the \code{ADCV}. Default is 15.
}
\item{unbiased}{
A logical value. If unbiased = TRUE, the unbiased estimator of squared
auto-distance covariance is returned.
Default value is FALSE.
}
}
\details{
Szekely et al. (2007) proposed distance covariance function between two random
vectors. Zhou (2012) extended this measure of dependence to a time series
framework by calling it auto-distance covariance function.

\code{\link{ADCV}} computes the sample auto-distance covariance function,
\eqn{V_X(\cdot)}, between \eqn{\{X_t\}} and \eqn{\{X_{t+j}\}}. Formal
definition of \eqn{V_X(\cdot)} can be found in Zhou (2012) and Fokianos and
Pitsillou (2017).

The empirical auto-distance covariance function, \eqn{\hat{V}_X(\cdot)},
is the non-negative square root defined by
\deqn{ \hat{V}_X^2(j) = \frac{1}{(n-j)^2}\sum_{r,l=1+j}^{n}{A_{rl}B_{rl}}, \quad 0 \leq j \leq (n-1)}{%
1/(n-j)^2 \sum_{r,l=1+j}^{n} A_{rl} B_{rl}, 0 \leq j \leq (n-1)}

and \eqn{\hat{V}_X^2(j) = \hat{V}_X^2(-j)}, for \eqn{-(n-1) \leq j < 0},
where \eqn{A=A_{rl}} and \eqn{B=B_{rl}} are Euclidean distances with elements
given by
\deqn{ A_{rl} = a_{rl} - \bar{a}_{r.} - \bar{a}_{.l} + \bar{a}_{..}
}
with \eqn{a_{rl}=|X_r-X_l|}, \eqn{\bar{a}_{r.}=\Bigl(\sum_{l=1+j}^{n}{a_{rl}}\Bigr)/(n-j)}{%
\bar{a}_{r.} = (\sum_{1+j}^{n} a_{rl})/(n-j)},
\eqn{\bar{a}_{.l}=\Bigl(\sum_{r=1+j}^{n}{a_{rl}}\Bigr)/(n-j)}{%
\bar{a}_{.l} = (\sum_{r=1+j}^{n} a_{rl})/(n-j)}
, \eqn{\bar{a}_{..}=\Bigl(\sum_{r,l=1+j}^{n}{a_{rl}}\Bigr)/(n-j)^2}{%
\bar{a}_{..}=(\sum_{r,l=1+j}^{n} a_{rl})/(n-j)^2
}.
\eqn{B_{rl}} is given analogously based on \eqn{b_{rl}=|Y_r-Y_l|}, where
\eqn{Y_t=X_{t+j}}. \eqn{X_t} and
\eqn{X_{t+j}} are independent if and only if \eqn{V_X^2(j)=0}.
See Fokianos and Pitsillou (2017) for more information on theoretical
properties of \eqn{V_X^2(\cdot)} including consistency.

If unbiased = TRUE, \code{\link{ADCV}} returns the unbiased estimator of
squared auto-distance covariance function,
\eqn{\tilde{V}_X^2(j)}, proposed by Szekely and Rizzo (2014).
In the context of time series data, this is given by
\deqn{  \tilde{V}_X^2(j) = \frac{1}{(n-j)(n-j-3)}\sum_{r\neq l}{\tilde{A}_{rl}\tilde{B}_{rl}},
}{%
\tilde{V}_X^2(j) = 1/((n-j)(n-j-3)) \sum_{r \neq l} \tilde{A}_{rl} \tilde{B}_{rl},
}
for \eqn{n > 3}, where \eqn{\tilde{A}_{rl}} is the \eqn{(r,l)}
element of the so-called U-centered matrix
\eqn{\tilde{A}}, defined by \deqn{ \tilde{A}_{rl} = \frac{1}{n-j-2}\sum_{t=1+j}^{n}{a_{rt}}-
\frac{1}{n-j-2}\sum_{s=1+j}^{n}{a_{sl}+\frac{1}{(n-j-1)(n-j-2)}\sum_{t,s=1+j}^{n}{a_{ts}}}, \quad i \neq j,
}{%
\tilde{A}_{rl} = 1/(n-j-2) \sum_{t=1+j}^{n} a_{rt}- 1/(n-j-2) \sum_{s=1+j}^{n} a_{sl}+ 1/( (n-j-1)(n-j-2) )
\sum_{t,s=1+j}^{n} a_{ts}, i \neq j,
}
with zero diagonal.

\code{\link{mADCV}} gives the auto-distance covariance function of a
multivariate time series.
}
\value{
A vector whose length is determined by \code{MaxLag} and contains the biased
estimator of \code{ADCV} or the unbiased estimator of squared \code{ADCV}.
}
\references{
Dominic, E, K. Fokianos and M. Pitsillou Maria (2019). An Updated Literature
Review of Distance Correlation and Its Applications to Time Series.
\emph{International Statistical Review}, 87, 237-262. .

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise dependence
in time series. \emph{Technometrics}, 159(2), 262-3270.

Huo, X. and G. J. Szekely. (2016). Fast Computing for Distance Covariance.
\emph{Technometrics}, 58, 435-447.

Pitsillou M. and Fokianos K. (2016). dCovTS: Distance Covariance/Correlation
for Time Series. \emph{R Journal}, 8, 324-340.

Szekely, G. J. and M. L. Rizzo (2014). Partial distance correlation with
methods for dissimilarities. \emph{The Annals of Statistics} 42, 2382-2412.

Szekely, G. J., M. L. Rizzo and N. K. Bakirov (2007). Measuring and testing
dependence by correlation of distances. \emph{The Annals of Statistics} 35, 2769-2794.

Zhou, Z. (2012). Measuring nonlinear dependence in time series,
a distance correlation approach.
\emph{Journal of Time Series Analysis} 33, 438-457.
}
\author{
Maria Pitsillou, Michail Tsagris and Konstantinos Fokianos.
}
\note{
Based on the definition of \eqn{\hat{V}_X(\cdot)}, we observe that
\eqn{\hat{V}^2_X(j)=\hat{V}^2_X(-j)}, and thus results based on negative
lags are omitted.
}


\seealso{
\code{\link{ADCF}}, \code{\link{mADCV}}
}

\examples{
x <- rnorm(500)
ADCV(x, 18)

ADCV(BJsales, 25)
}
\keyword{ ts }
