\name{mADCVtest}
\alias{mADCVtest}
\title{
Distance covariance test of independence in multivariate time series
}
\description{
A test of independence based on auto-distance covariance matrix in multivariate time series proposed by Fokianos a
nd Pitsillou (2017).
}
\usage{
mADCVtest(x, type = c("truncated", "bartlett", "daniell", "QS", "parzen"), p,
          b = 0, parallel = FALSE, bootMethod = c("Wild Bootstrap",
          "Independent Bootstrap"))
}
\arguments{
\item{x}{
Multivariate time series.
}
\item{type}{
A character string which indicates the smoothing kernel. Possible choices are 'truncated' (the default), 'bartlett',
'daniell', 'QS', 'parzen'.
}
\item{p}{
The bandwidth, whose choice is determined by \eqn{p=cn^{\lambda}} for \eqn{c > 0} and \eqn{\lambda \in (0,1)}.
}
\item{b}{
The number of bootstrap replicates of the test statistic. It is a positive integer. If b=0 (the default), then no
p-value is returned.
}
\item{parallel}{
A logical value. By default, parallel=FALSE. If parallel=TRUE, bootstrap computation is distributed to multiple cores,
which typically is the maximum number of available CPUs and is detecting directly from the function.
}
\item{bootMethod}{
A character string indicating the method to use for obtaining the empirical p-value of the test.
Possible choices are "Wild Bootstrap" (the default) and "Independent Bootstrap".
}
}
\details{
\code{\link{mADCVtest}} tests whether the vector series are independent and identically distributed (i.i.d). The p-value
of the test is obtained via resampling scheme. Possible choices are the independent wild bootstrap (Dehling and Mikosch, 1994; Shao, 2010;
Leucht and Neumann, 2013) and independent bootstrap, with \code{b} replicates. The observed statistic is
\deqn{
\sum_{j=1}^{n-1}(n-j)k^2(j/p)\mbox{tr}\{\hat{V}^{*}(j)\hat{V}(j)\}
}{%
\sum_{j=1}^{n-1} (n-j)k^2(j/p) tr {\hat{V}^{*}(j)\hat{V}(j) }
}
where \eqn{\hat{V}^{*}(\cdot)} denotes the complex conjugate matrix of \eqn{\hat{V}(\cdot)} obtained from \code{\link{mADCV}}, and
\eqn{\mbox{tr}\{A\}} denotes the trace of a matrix \eqn{A}, which is the sum of the diagonal elements of \eqn{A}. \eqn{k(\cdot)}
is a kernel function computed by \code{\link{kernelFun}} and \code{p} is a bandwidth or lag order whose choice is further discussed
in Fokianos and Pitsillou (2017).

Under the null hypothesis of independence and some further assumptions about the kernel function \eqn{k(\cdot)}, the standardized
version of the test statistic follows \eqn{N(0,1)} asymptotically and it is consistent. More details of the asymptotic properties
of the statistic can be found in Fokianos and Pitsillou (2017).

\code{\link{mADCFtest}} performs the same test based on the distance correlation matrix \code{\link{mADCF}}.
}
\value{
An object of class \code{htest} which is a list including:
\item{method}{
The description of the test.
}
\item{statistic}{
The observed value of the test statistic.
}
\item{replicates}{
Bootstrap replicates of the test statistic (if \eqn{b=0} then \code{replicates}=NULL).
}
\item{p.value}{
The p-value of the test (if \eqn{b=0} then \code{p.value}=NA).
}
\item{bootMethod}{
The method followed for computing the p-value of the test.
}
\item{data.name}{
The description of the data (data name, kernel type, \code{type}, bandwidth, \code{p}, and the number of
bootstrap replicates \code{b}).
}
}
\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature Review of Distance Correlation and Its
Applications to Time Series. \emph{International Statistical Review}, 87, 237-262.

Dehling, H. and T. Mikosch (1994). Random quadratic forms and the bootstrap for U-statistics.
\emph{Journal of Multivariate Analysis}, 51, 392-413.

Fokianos K. and Pitsillou M. (2018). Testing independence for multivariate time series via the auto-distance correlation matrix.
\emph{Biometrika}, 105, 337-352.

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise dependence in time series.
\emph{Technometrics}, 159, 262-3270.

Huo, X. and G. J. Szekely. (2016). Fast Computing for Distance Covariance.
\emph{Technometrics}, 58, 435-447.

Leucht, A. and M. H. Neumann (2013). Dependent wild bootstrap for degenerate U- and V- statistics.
\emph{Journal of Multivariate Analysis}, 117, 257-280.

Pitsillou M. and Fokianos K. (2016). dCovTS: Distance Covariance/Correlation for Time Series.
\emph{R Journal}, 8, 324-340.

Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105, 218-235.
}
\author{
Maria Pitsillou, Michail Tsagris and Konstantinos Fokianos.
}
\note{
The computation of the test statistic is only based on the biased estimator of auto-distance covariance matrix.
}

\seealso{
\code{\link{mADCV}}, \code{\link{mADCF}}, \code{\link{mADCFtest}}
}
\examples{
x <- matrix( rnorm(200), ncol = 2 )
n <- length(x)
c <- 3
lambda <- 0.1
p <- ceiling(c * n^lambda)
mF <- mADCVtest(x, type = "bar", p = p, b = 500, parallel = FALSE)
}
