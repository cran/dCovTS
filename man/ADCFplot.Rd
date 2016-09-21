\name{ADCFplot}
\alias{ADCFplot}

\title{
Auto-distance correlation plot
}
\description{
The function plots the estimated auto-distance correlation function obtained by \code{\link{ADCF}}.
}
\usage{
ADCFplot(x, MaxLag = 15, ylim = NULL, main = NULL, bootMethod = c("Wild Bootstrap", 
         "Subsampling", "Independent Bootstrap"), b = 499)
}

\arguments{
  \item{x}{
numeric vector or univariate time series.
}
  \item{MaxLag}{
maximum lag order at which to plot \code{\link{ADCF}}. Default is 15.
}
  \item{ylim}{
numeric vector of length 2 indicating the \code{y} limits of the plot. The default value, NULL, indicates that
the range \eqn{(0,v)}, where \eqn{v} is the maximum number between 1 and the empirical critical values, should be used.
}
  \item{main}{
title of the plot.
}
  \item{bootMethod}{
character string indicating the method to use for obtaining the 95\% critical values. Possible choices are "Wild Bootstrap" (the default), "Independent Bootstrap"
and "Subsampling".
}
  \item{b}{
the number of bootstrap replications for constructing the 95\% empirical critical values. Default is 499.
}
}
\details{
Fokianos and Pitsillou (2016) showed that the sample auto-distance covariance function \code{\link{ADCV}} (and thus \code{\link{ADCF}}) can be expressed as a V-statistic of order two, which
under the null hypothesis of independence is degenerate. Thus, constructing a plot analogous to the traditional autocorrelation
plot where the confidence intervals are obtained simultaneously, turns to be a complicated task. To overcome this issue, the 95\% 
confidence intervals shown in the plot (dotted blue horizontal line) are computed simultaneously via Monte Carlo simulation, and in
particular via the independent wild bootstrap approach (Dehling and Mikosch, 1994; Shao, 2010; Leucht and Neumann, 2013). The reader is referred to Fokianos and Pitsillou (2016, Section 6.2) 
for the steps followed. \code{\link{mADCFplot}} returns an analogous plot of the estimated auto-distance correlation function for a 
multivariate time series.


One can also compute the pairwise 95\% critical values via the subsampling approach suggested by Zhou (2012, Section 5.1).
That is, the critical values are obtained at each lag separately. The block size of the procedure is based on the minimum volatility
method proposed by Politis et al. (1999, Section 9.4.2).  In addition, the function provides the ordinary independent bootstrap methodology to derive
simultaneous 95\% critical values.
}
\value{
A plot of the estimated \code{\link{ADCF}} values. It also returns a list with
\item{ADCF}{
The sample auto-distance correlation function for all lags specified by \code{MaxLag}.
}
\item{bootMethod}{
The method followed for computing the 95\% confidence intervals of the plot.
}
\item{critical.value}{
The critical value shown in the plot.
}
}
\references{
Dehling, H. and T. Mikosch (1994). Random quadratic forms and the bootstrap for U-statistics. 
\emph{Journal of Multivariate Analysis} \eqn{\textbf{51}}{51}, 392-413,
\url{http://dx.doi.org/10.1006/jmva.1994.1069}

Fokianos K. and M. Pitsillou (2016). Testing pairwise independence for multivariate
time series by the auto-distance correlation matrix. Submitted for publication.

Leucht, A. and M. H. Neumann (2013). Dependent wild bootstrap for degenerate U- and V- statistics.
\emph{Journal of Multivariate Analysis} \eqn{\textbf{117}}{117}, 257-280,
\url{http://dx.doi.org/10.1016/j.jmva.2013.03.003}.

Politis, N. P., J. P. Romano and M. Wolf (1999). \emph{Subsampling}. New York: Springer.

Shao, X. (2010). The dependent wild bootstrap.
\emph{Journal of the American Statistical Association} \eqn{\textbf{105}}{105}, 218-235,
\url{http://dx.doi.org/10.1198/jasa.2009.tm08744}.

Zhou, Z. (2012). Measuring nonlinear dependence in time series, a distance correlation approach.
\emph{Journal of Time Series Analysis} \eqn{\textbf{33}}{33}, 438-457,
\url{http://dx.doi.org/10.1111/j.1467-9892.2011.00780.x}.
}
\note{
When the critical values are obtained via the Subsampling methodology, the function returns a plot that starts from lag 1.

The function plots only the biased estimator of ADCF.
}

\seealso{
\code{\link{ADCF}}, \code{\link{ADCV}}, \code{\link{mADCFplot}}
}
\examples{
\dontrun{ADCFplot(rnorm(100),ylim=c(0,0.4),bootMethod="Subs")}

ADCFplot(mdeaths,bootMethod="Wild",b=100)

ADCFplot(mdeaths,bootMethod="Indep",b=100)
}
\keyword{ plot }
