\name{UnivTest}
\alias{UnivTest}

\title{
Testing for independence in univariate time series
}
\description{
A test of pairwise independence for univariate time series.
}
\usage{
UnivTest(x, type = c("truncated", "bartlett", "daniell", "QS", "parzen"), 
         testType = c("covariance", "correlation"), p, b = 0, parallel = FALSE, 
         bootMethod = c("Wild Bootstrap", "Independent Bootstrap"))
}

\arguments{
  \item{x}{
numeric vector or univariate time series.
}
  \item{type}{
character string which indicates the smoothing kernel. Possible choices are 'truncated' (the default), 'bartlett', 'daniell', 'QS', 'parzen'.
}
  \item{testType}{
character string indicating the type of the test to be used. Allowed values are 'covariance' (default) for using
the distance covariance function and 'correlation' for using the distance correlation function.
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
  \item{bootMethod}{
character string indicating the method to use for obtaining the empirical p-value of the test. 
Possible choices are "Wild Bootstrap" (the default) and "Independent Bootstrap"
}
}
\details{
\code{UnivTest} performs a test on the null hypothesis of independence in univariate time series.
The p-value of the test is obtained via resampling method. Possible choices are the independent wild bootstrap (Dehling and Mikosch, 1994; 
Shao, 2010; Leucht and Neumann, 2013) (default option) and the ordinary independent bootstrap,
with \code{b} replicates. If \code{typeTest} = 'covariance' then, the observed statistic is
\deqn{
\sum_{j=1}^{n-1}{(n-j)k^2(j/p)\hat{V}^2_X(j)},
}{%
\sum_{j=1}^{n-1} (n-j)k^2(j/p)\hat{V}^2_X(j) ,
}
otherwise
\deqn{
\sum_{j=1}^{n-1}{(n-j)k^2(j/p)\hat{R}^2_X(j)},
}{%
\sum_{j=1}^{n-1} (n-j)k^2(j/p)\hat{R}^2_X(j) ,
}
where \eqn{k(\cdot)} is a kernel function computed by \code{\link{kernelFun}} and \code{p} is a bandwidth or lag order whose
choice is further discussed in Fokianos and Pitsillou (2016).

Under the null hypothesis of independence and some further assumptions about the kernel function \eqn{k(\cdot)}, the standardized version of the
test statistic follows \eqn{N(0,1)} asymptotically and it is consistent. More 
details of the asymptotic properties of the statistic can be found in Fokianos and Pitsillou (2016). 
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
  \item{bootMethod}{
The method followed for computing the p-value of the test.
}
  \item{data.name}{
description of data (the data name, kernel type, \code{type}, bandwidth, \code{p}, and the number of bootstrap replicates \code{b}).
}

}
\references{
Dehling, H. and T. Mikosch (1994). Random quadratic forms and the bootstrap for U-statistics. 
\emph{Journal of Multivariate Analysis} \eqn{\textbf{51}}{51}, 392-413,
\url{http://dx.doi.org/10.1006/jmva.1994.1069}

Fokianos K. and M. Pitsillou (2016). Consistent testing for pairwise dependence in time series. \emph{Technometrics},
\url{http://dx.doi.org/10.1080/00401706.2016.1156024}.

Leucht, A. and M. H. Neumann (2013). Dependent wild bootstrap for degenerate U- and V- statistics.
\emph{Journal of Multivariate Analysis} \eqn{\textbf{117}}{117}, 257-280,
\url{http://dx.doi.org/10.1016/j.jmva.2013.03.003}.

Shao, X. (2010). The dependent wild bootstrap.
\emph{Journal of the American Statistical Association} \eqn{\textbf{105}}{105}, 218-235,
\url{http://dx.doi.org/10.1198/jasa.2009.tm08744}.
}

\author{
Maria Pitsillou and Konstantinos Fokianos
}
\note{
The observed statistics of the tests are only based on the biased estimators of distance covariance and correlation functions.
}

\seealso{
\code{\link{ADCF}} \code{\link{ADCV}}
}
\examples{
x<-rnorm(100)
n <- length(x)
c <- 1
lambda <- 1/5
p <- ceiling(c*n^lambda)
\dontrun{
mW=UnivTest(x,type="bar",testType="covariance",p=p,b=499,parallel=TRUE,bootMethod="Wild")
mI=UnivTest(x,type="bar",testType="covariance",p=p,b=499,parallel=TRUE,bootMethod="Indep")
}

data <- tail(ibmSp500[,2],100)
n2 <- length(data)
c2 <- 3
lambda2 <- 0.1
p2 <- ceiling(c2*n2^lambda2)
\dontrun{
testCov=UnivTest(data,type="par",testType="covariance",p=p2,b=499,parallel=TRUE)
testCor=UnivTest(data,type="par",testType="correlation",p=p2,b=499,parallel=TRUE)
}
}

\keyword{ htest }
