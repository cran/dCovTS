\name{mADCFplot}
\alias{mADCFplot}

\title{
Distance cross-correlation plot
}
\description{
The function computes and plots the estimator of the auto-distance correlation
matrix \code{\link{mADCF}}.
}
\usage{
mADCFplot(x, MaxLag = 15, alpha = 0.05, b = 499,
          bootMethod = c("Wild Bootstrap", "Independent Bootstrap"),
          ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A multivariate time series.
}
\item{MaxLag}{
The maximum lag order at which to plot \code{\link{mADCF}}. Default is 15.
}
\item{alpha}{
The significance level used to construct the \eqn{(1-\alpha)}\% empirical
critical values.
}
\item{b}{
The number of bootstrap replications for constructing the \eqn{(1-\alpha)}\%
empirical critical values. Default is 499.
}
\item{bootMethod}{
A character string indicating the method to use for obtaining the
\eqn{(1-\alpha)}\% critical values. Possible choices are
"Wild Bootstrap" (the default) and "Independent Bootstrap".
}
\item{ylim}{
A numeric vector of length 2 indicating the \code{y} limits of the plot.
The default value, NULL, indicates that the range \eqn{(0,v)}, where
\eqn{v} is the maximum number between 1 and the empirical critical values,
should be used.
}
}
\details{
The \eqn{(1-\alpha)}\% confidence intervals shown in the plot
(dotted blue horizontal line) are computed simultaneously based
on the independent wild bootstrap approach (Dehling and Mikosch,
1994; Shao, 2010; Leucht and Neumann, 2013), since the
elements of \code{\link{mADCV}} (and thus \code{\link{mADCF}}) can be
expressed as degenerate V-statistics of order 2.
More details can be found in Fokianos and Pitsillou (2017).

In addition, \code{\link{mADCFplot}} provides the option of independent
bootstrap to compute the simultaneous
\eqn{(1-\alpha)}\% critical values.
}
\value{
A plot of the estimated \code{\link{mADCF}} matrices. The function also
returns a list including
\item{matrices}{
Sample distance correlation matrices starting from lag 0.
}
\item{bootMethod}{
The method followed for computing the \eqn{(1-\alpha)}\% confidence
intervals of the plot.
}
\item{critical.value}{
The critical value shown in the plot.
}
}
\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature
Review of Distance Correlation and Its Applications to Time Series.
\emph{International Statistical Review}, 87, 237-262.

Dehling, H. and T. Mikosch (1994). Random quadratic forms and the
bootstrap for U-statistics. \emph{Journal of Multivariate Analysis}, 51, 392-413.

Fokianos K. and Pitsillou M. (2018). Testing independence for multivariate time
series via the auto-distance correlation matrix. \emph{Biometrika}, 105, 337-352.

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise dependence
in time series. \emph{Technometrics}, 159, 262-3270.

Huo, X. and G. J. Szekely. (2016). Fast Computing for Distance Covariance.
\emph{Technometrics}, 58, 435-447.

Leucht, A. and M. H. Neumann (2013). Dependent wild bootstrap for degenerate
U- and V- statistics. \emph{Journal of Multivariate Analysis}, 117, 257-280.

Pitsillou M. and Fokianos K. (2016). dCovTS: Distance Covariance/Correlation
for Time Series. \emph{R Journal}, 8, 324-340.

Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American
Statistical Association}, 105, 218-235.
}
\author{
Maria Pitsillou and Konstantinos Fokianos.
}
\note{
The function plots only the biased estimator of ADCF matrix.

%Due to restrictions in the use of parallel on CRAN, we removed the example.
%Two examples of how this function works are
%"x <- matrix( rnorm(200), ncol = 2 )";
%"mADCFplot(x, 12, ylim = c(0, 0.5) )";
%"mADCFplot(x, 12, b = 100)".
}


\seealso{
\code{\link{mADCF}}, \code{\link{mADCV}}
}
\examples{
\donttest{
### x <- matrix( rnorm(200), ncol = 2 )
### mADCFplot(x, 12, ylim = c(0, 0.5) )
### mADCFplot(x, 12, b = 100)
}
}
\keyword{ plot }
