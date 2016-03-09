\name{ADCV}
\alias{ADCV}

\title{
Auto-distance Covariance Function
}
\description{
Computes the auto-distance covariance function of a univariate time series.
}
\usage{
ADCV(x, MaxLag)
}
\arguments{
  \item{x}{
numeric vector or univariate time series.
}
  \item{MaxLag}{
maximum lag order at which to calculate the \code{ADCV}. Default is 15.
}
}
\details{
Szekely et al. (2007) recently proposed distance covariance function between two random vectors. Zhou (2012) extended this measure of dependence
to a time series framework by calling it auto-distance covariance function.
 
\code{ADCV} computes the sample auto-distance covariance function, \eqn{V_X(\cdot)}, between \eqn{\{X_t\}} and \eqn{\{X_{t-|j|}\}}.
Formal definition of \eqn{V_X(\cdot)} can be found in Zhou (2012) and Fokianos and Pitsillou (2016).

%The auto-distance covariance function, \eqn{V_X(j)}, between \eqn{\{X_t\}} and \eqn{\{X_{t-|j|}\}} is defined by
%\deqn{ V_X^2(j) = \int_{\Re^2} |\sigma_j(u,v)|^2 W(u,v)dudv
% }{ (non-Latex version) }
%for \eqn{j = 0, \pm 1, \pm 2, \dots} and \eqn{(u,v) \in \Re^2}. The dependence measure \eqn{\sigma_j(u,v)} was first defined by Hong(1999) and is given by
%\deqn{ \sigma_j(u,v)= Cov(e^{iuX_t},e^{iuX_{t-|j|}})
%}
%and the weight function (Szekely et al, 2007;Zhou, 2012) is given by
%\deqn{ W(u,v)=\frac{1}{\pi^2 |u|^2|v|^2}.
%}
The empirical auto-distance covariance function, \eqn{\hat{V}_X(\cdot)}, is the non-negative square root defined by
\deqn{ \hat{V}_X^2(j) = \frac{1}{(n-|j|)^2}\sum_{r,l=1+|j|}^{n}{A_{rl}B_{rl}}, ~~~~~ j=0, \pm 1, \pm 2, \dots.
}

where \eqn{A=A_{rl}} and \eqn{B=B_{rl}} are Euclidean distances with elements given by
\deqn{ A_{rl} = a_{rl} - \bar{a}_{r.} - \bar{a}_{.l} + \bar{a}_{..} 
}
with \eqn{a_{rl}=|X_r-X_l|}, \eqn{\bar{a}_{r.}=\Bigl(\sum_{l=1+|j|}^{n}{a_{rl}}\Bigr)/(n-|j|)},
\eqn{\bar{a}_{.l}=\Bigl(\sum_{r=1+|j|}^{n}{a_{rl}}\Bigr)/(n-|j|)}, \eqn{\bar{a}_{..}=\Bigl(\sum_{r,l=1+|j|}^{n}{a_{rl}}\Bigr)/(n-|j|)^2}.
\eqn{B_{rl}} is given analogously based on \eqn{b_{rl}=|Y_r-Y_l|}, where \eqn{Y_t=X_{t-|j|}}. \eqn{X_t} and \eqn{X_{t-|j|}} are independent if and only if \eqn{V_X^2(j)=0}. See Fokianos and Pitsillou (2016) for more information on theoretical properties
of \eqn{V_X^2(\cdot)} including consistency. 
 
\code{\link{mADCV}} gives the auto-distance covariance function of a multivariate time series.
}
\value{
Returns a vector, whose length is determined by \code{MaxLag}, and contains the estimated \code{ADCV}. 
}
\references{
Fokianos K. and M. Pitsillou (2016). Consistent testing for pairwise dependence in time series. \emph{Technometrics},
\url{http://dx.doi.org/10.1080/00401706.2016.1156024}.
 
Szekely, G. J., M. L. Rizzo and N. K. Bakirov (2007). Measuring and testing dependence 
by correlation of distances. \emph{The Annals of Statistics} \eqn{\textbf{35}}, 2769-2794,
\url{http://dx.doi.org/10.1214/009053607000000505}.

Zhou, Z. (2012). Measuring nonlinear dependence in time series, a distance correlation approach.
\emph{Journal of Time Series Analysis} \eqn{\textbf{33}}, 438-457,
\url{http://dx.doi.org/10.1111/j.1467-9892.2011.00780.x}.

}
\author{
Maria Pitsillou and Konstantinos Fokianos
}
\note{
Based on the definition of \eqn{V_X(\cdot)}, we observe that \eqn{V^2_X(j)=V^2_X(-j)}, and thus results based on negative lags are omitted.

\code{ADCV} returns \eqn{\hat{V}_X(\cdot)} and not \eqn{\hat{V}^2_X(\cdot)} \eqn{(\code{ADCV}^2)}.
}

\seealso{
\code{\link{ADCF}}, \code{\link{mADCV}}
}

\examples{
x <- rnorm(500)
ADCV(x,18) 

ADCV(BJsales,25)
}
\keyword{ ts }

