\name{mADCFplot}
\alias{mADCFplot}

\title{
Distance cross-correlation plot
}
\description{
The function computes and plots the estimator of the distance correlation matrix \code{\link{mADCF}}.
}
\usage{
mADCFplot(x, MaxLag = 15, ylim= NULL, b = 499)
}

\arguments{
  \item{x}{
multivariate time series.
}
  \item{MaxLag}{
maximum lag order at which to plot \code{\link{mADCF}}. Default is 15.
}
  \item{ylim}{
numeric vector of length 2 indicating the \code{y} limits of the plot. The default value, NULL, indicates that
the range \eqn{(0,v)}, where \eqn{v} is the maximum number between 1 and the empirical critical values, should be used.
}
  \item{b}{
the number of Wild bootstrap replications for constructing the 95\% empirical confidence intervals. Default is 499.
}
}
\details{
The 95\% confidence intervals shown in the plot (dotted blue horizontal line) are computed simultaneously based on the Independent Wild Bootstrap approach (Shao, 2010; Leucht and Neumann, 2013),
since the elements of \code{\link{mADCV}} (and thus \code{\link{mADCF}}) can be expressed as degenerate V-statistics of order 2. 
More details can be found in Fokianos and Pitsillou (2016).
}
\value{
A plot of the estimated \code{\link{mADCF}} matrices. The function also returns a list with
\item{matrices}{
Sample distance correlation matrices starting from lag 0.
}
\item{critical.value}{
The critical value shown in the plot.
}
}
\references{
Fokianos K. and M. Pitsillou (2016). On multivariate auto-distance covariance and correlation functions. 
Submitted for publication.

Leucht, A. and M. H. Neumann (2013). Dependent wild bootstrap for degenerate U- and V- statistics.
\emph{Journal of Multivariate Analysis} \eqn{\textbf{117}}, 257-280,
\url{http://dx.doi.org/10.1016/j.jmva.2013.03.003}.

Shao, X. (2010). The dependent wild bootstrap.
\emph{Journal of the American Statistical Association} \eqn{\textbf{105}}, 218-235,
\url{http://dx.doi.org/10.1198/jasa.2009.tm08744}.
}
\author{
Maria Pitsillou and Konstantinos Fokianos
}


\seealso{
\code{\link{mADCF}}, \code{\link{mADCV}}
}
\examples{
x <- MASS::mvrnorm(100,rep(0,3),diag(3))
\dontrun{mADCFplot(x,18,ylim=c(0,0.5))}
y <- MASS::mvrnorm(100,rep(0,6),diag(6))
\dontrun{mADCFplot(y,b=100)}

deaths <- cbind(mdeaths,fdeaths)
\dontrun{mADCFplot(deaths)}
}
\keyword{ plot }
