\name{mADCVtest}
\alias{mADCVtest}
\title{
Distance covariance test of independence in multivariate time series
}
\description{
A test of independence based on distance covariance function in multivariate time series proposed by Fokianos and Pitsillou (2016).
}
\usage{
mADCVtest(x, type = c("truncated", "bartlett", "daniell", "QS", "parzen"), p,
          b = 0, parallel = FALSE)
}

\arguments{
  \item{x}{
multivariate time series.
}
  \item{type}{
character string which indicates the smoothing kernel. Possible choices are 'truncated' (the default), 'bartlett', 'daniell', 'QS', 'parzen'.
}
  \item{p}{
bandwidth, whose choice is determined by \eqn{p=cn^{\lambda}} for \eqn{c > 0} and \eqn{\lambda \in (0,1)}.
}
  \item{b}{
the number of bootstrap replicates of the test statistic. It is a positive integer. If b=0 (the default), then no p-value is returned.
}
  \item{parallel}{
logical value. By default, parallel=FALSE. If parallel=TRUE, bootstrap computation is distributed to multiple cores, which 
typically is the maximum number of available CPUs and is detecting directly from the function.
}
}
\details{
\code{mADCVtest} performs a test of multivariate independence. In particular, the function tests whether there are no auto- and cross-correlations
in a multivariate time series. The p-value of the test is obtained via the independent wild bootstrap (Shao, 2010; Leucht and Neumann, 2013),
with \code{b} replicates. The observed statistic is
\deqn{
\sum_{j=1}^{n-1}(n-j)k^2(j/p)\mbox{tr}\{\hat{V}^{*}(j)\hat{V}(j)\}
}
where \eqn{\hat{V}^{*}(\cdot)} denotes the complex conjugate matrix of \eqn{\hat{V}(\cdot)} obtained from \code{\link{mADCV}}, and
\eqn{\mbox{tr}\{A\}} denotes the trace of a matrix \eqn{A}, which is the sum of the diagonal elements of \eqn{A}. \eqn{k(\cdot)} is a kernel function computed by \code{\link{kernelFun}} and \code{p} is a bandwidth or lag order whose
choice is further discussed in Fokianos and Pitsillou (2016).  

Under the null hypothesis of independence and some further assumptions about the kernel function \eqn{k(\cdot)}, the standardized version of the
test statistic follows \eqn{N(0,1)} asymptotically and it is consistent. More 
details of the asymptotic properties of the statistic can be found in Fokianos and Pitsillou (2016). 

\code{\link{mADCFtest}} performs the same test based on the distance correlation matrix \code{\link{mADCF}}. 
}

\value{
An object of class \code{htest} which is a list containing:
  \item{method}{
description of test.
}
  \item{statistic}{
the observed value of the test statistic.
}
  \item{replicates}{
bootstrap replicates of the test statistic (if \eqn{b=0} then \code{replicates}=NULL).
}
  \item{p.value}{
p-value of the test (if \eqn{b=0} then \code{p.value}=NA).
}
  \item{data.name}{
description of data (data name, kernel type, \code{type}, bandwidth, \code{p}, and the number of bootstrap replicates \code{b}).
}
}

\references{
Fokianos K. and M. Pitsillou (2016). On multivariate auto-distance covariance and correlation functions. Submitted for publication.

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
\code{\link{mADCV}}, \code{\link{mADCF}}, \code{\link{mADCFtest}}
}
\examples{
x<-MASS::mvrnorm(500,rep(0,2),diag(2))
n <- length(x)
c <- 3
lambda <- 0.1
p <- ceiling(c*n^lambda)
\dontrun{
mT=mADCVtest(x,type="bar",p=p,b=499,parallel=TRUE)
mF=mADCVtest(x,type="bar",p=p,b=499,parallel=FALSE)}
}

\keyword{ htest }
