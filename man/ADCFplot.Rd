\name{ADCFplot}
\alias{ADCFplot}
\title{
Auto-distance correlation plot
}
\description{
The function plots the estimated auto-distance correlation function obtained by \code{\link{ADCF}}
and provides confindence intervals by employing three bootstrap based methods.
}
\usage{
ADCFplot(x, MaxLag = 15, alpha = 0.05, b = 499, bootMethod =
c("Wild Bootstrap", "Subsampling", "Independent Bootstrap"), ylim = NULL, main = NULL)
}

\arguments{
\item{x}{
A numeric vector or univariate time series.
}
\item{MaxLag}{
The maximum lag order at which to plot \code{\link{ADCF}}. Default is 15.
}
\item{alpha}{
The significance level used to construct the \eqn{(1-\alpha)}\% empirical critical values.
}
\item{b}{
The number of bootstrap replications for constructing the \eqn{(1-\alpha)}\% empirical critical values.
Default is 499.
}
\item{bootMethod}{
A character string indicating the method to use for obtaining the \eqn{(1-\alpha)}\% critical values. Possible
choices are "Wild Bootstrap" (the default), "Independent Bootstrap" and "Subsampling".
}
\item{ylim}{
A numeric vector of length 2 indicating the \code{y} limits of the plot. The default value, NULL, indicates
that the range \eqn{(0,v)}, where \eqn{v} is the maximum number between 1 and the empirical critical values,
should be used.
}
\item{main}{
The title of the plot.
}
}
\details{
Fokianos and Pitsillou (2018) showed that the sample auto-distance covariance function \code{\link{ADCV}}
(and thus \code{\link{ADCF}}) can be expressed as a V-statistic of order two, which under the null hypothesis
of independence is degenerate. Thus, constructing a plot analogous to the traditional autocorrelation plot where
the confidence intervals are obtained simultaneously, turns to be a complicated task. To overcome this issue, the
\eqn{(1-\alpha)}\% confidence intervals shown in the plot (dotted blue horizontal line) are computed simultaneously
via Monte Carlo simulation, and in particular via the independent wild bootstrap approach
(Dehling and Mikosch, 1994; Shao, 2010; Leucht and Neumann, 2013). The reader is referred to Fokianos and
Pitsillou (2018) for the steps followed. \code{\link{mADCFplot}} returns an analogous plot of the estimated
auto-distance correlation function for a multivariate time series.

One can also compute the pairwise \eqn{(1-\alpha)}\% critical values via the subsampling approach suggested by
Zhou (2012, Section 5.1).That is, the critical values are obtained at each lag separately. The block size of the
procedure is based on the minimum volatility method proposed by Politis et al. (1999, Section 9.4.2). In addition,
the function provides the ordinary independent bootstrap methodology to derive simultaneous \eqn{(1-\alpha)}\% critical values.
}
\value{
A plot of the estimated \code{\link{ADCF}} values. It also returns a list including:
\item{ADCF}{
The sample auto-distance correlation function for all lags specified by \code{MaxLag}.
}
\item{bootMethod}{
The method followed for computing the \eqn{(1-alpha)}\% confidence intervals of the plot.
}
\item{critical.value}{
The critical value shown in the plot.
}
}
\references{
Dehling, H. and T. Mikosch (1994). Random quadratic forms and the bootstrap for U-statistics.
\emph{Journal of Multivariate Analysis}, 51, 392-413.

Dominic, E, K. Fokianos and M. Pitsillou Maria (2019). An Updated Literature Review of Distance Correlation
and Its Applications to Time Series. \emph{International Statistical Review}, 87, 237-262.

Fokianos K. and Pitsillou M. (2018). Testing independence for multivariate time series via the auto-distance
correlation matrix. \emph{Biometrika}, 105, 337-352.

Leucht, A. and M. H. Neumann (2013). Dependent wild bootstrap for degenerate U- and V- statistics.
\emph{Journal of Multivariate Analysis}, 117, 257-280.

Pitsillou M. and Fokianos K. (2016). dCovTS: Distance Covariance/Correlation for Time Series.
R Journal, 8, 324-340.

Politis, N. P., J. P. Romano and M. Wolf (1999). \emph{Subsampling}. New York: Springer.

Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105 218-235.

Zhou, Z. (2012). Measuring nonlinear dependence in time series, a distance correlation approach.
\emph{Journal of Time Series Analysis}, 33, 438-457.
}
\author{
Maria Pitsillou, Michail Tsagris and Konstantinos Fokianos.
}
\note{
When the critical values are obtained via the Subsampling methodology, the function returns a plot that
starts from lag 1.

The function plots only the biased estimator of ADCF.
}

\seealso{
\code{\link{ADCF}}, \code{\link{ADCV}}, \code{\link{mADCFplot}}
}
\examples{
\dontrun{
ADCFplot(rnorm(100), bootMethod = "Subs")
}
}
\keyword{ plot }
