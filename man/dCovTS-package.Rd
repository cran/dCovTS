\name{dCovTS-package}
\alias{dCovTS-package}
\alias{dCovTS}
\docType{package}
\title{
Distance Covariance and Correlation Theory for Time Series
}
\description{
Computing and plotting the distance covariance and correlation
function of a univariate or a multivariate time series.
Both versions of biased and unbiased estimators of distance covariance
and correlation are provided. Test statistics for testing pairwise independence
are also implemented. Some data sets are also included.
}
\details{
\tabular{ll}{
Package: \tab dCovTS\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2022-06-27\cr
License: \tab GPL(>=2)\cr
}
}
\author{
Michail Tsagris, Maria Pitsillou and Konstantinos Fokianos.
}
\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature
Review of Distance Correlation and Its Applications to Time Series.
\emph{International Statistical Review}, 87, 237-262.

Fokianos, K. and M. Pitsillou. (2018). Testing independence for multivariate
time series via the auto-distance correlation matrix.
\emph{Biometrika}, 105, 337-352.

Fokianos, K. and M. Pitsillou. (2017). Consistent testing for pairwise
dependence in time series.
\emph{Technometrics}, 159, 262-3270.

Dehling, H. and T. Mikosch. (1994). Random quadratic forms and the bootstrap
for U-statistics.
\emph{Journal of Multivariate Analysis} 51, 392-413.

Hong, Y. (1999). Hypothesis testing in time series via the empirical
characteristic function: A generalized spectral density approach.
\emph{Journal of the American Statistical Association}, 94, 1201-1220.

Hong, Y. (1996). Consistent testing for serial correlation of unknown form.
\emph{Econometrica}, 64, 837-864.

Huo, X. and G. J. Szekely. (2016). Fast Computing for Distance Covariance.
\emph{Technometrics}, 58, 435-447.

Leucht, A. and M. H. Neumann (2013). Dependent wild bootstrap for degenerate
U- and V- statistics.
\emph{Journal of Multivariate Analysis}, 117, 257-280.

Pitsillou, M. and K. Fokianos. (2016). dCovTS: Distance Covariance/Correlation
for Time Series.
\emph{R Journal}, 8, 324-340.

Politis, N. P., J. P. Romano and M. Wolf (1999).
\emph{Subsampling}. New York: Springer.

Shao, X. (2010). The dependent wild bootstrap.
\emph{Journal of the American Statistical Association}, 105, 218-235.

Shumway, R. H. and D. S. Stoffer (2011).
\emph{Time Series Analysis and Its Applications With R Examples}.
New York: Springer. Third Edition.

Szekely, G. J. and M. L. Rizzo (2014). Partial distance correlation with
methods for dissimilarities. \emph{The Annals of Statistics}, 42, 2382-2412.

Szekely, G. J., M. L. Rizzo and N. K. Bakirov (2007). Measuring and testing
dependence by correlation of distances.
\emph{The Annals of Statistics}, 35, 2769-2794, .

Tsay, R. S. (2014). \emph{Multivariate Time Series Analysis with R and
Financial Applications}. Hoboken, NJ: Wiley.

Tsay, R. S. (2010). \emph{Analysis of Financial Time Series}.
Hoboken, NJ: Wiley. Third edition.

Zhou, Z. (2012). Measuring nonlinear dependence in time series, a distance
correlation approach.
\emph{Journal of Time Series Analysis}, 33, 438-457.
}

\note{
\bold{Disclaimer:} Dr Maria Pitsillou is the actual creator of this package.
Dr Tsagris is the current maintainer.

Improvements: We have modified the codes to run faster, we included the
packages Rfast and Rfast2 for fast computations and the "dcov" package that
allows for extremely fast computations of the distance correlation/covariance
with univariate data.
}

\keyword{ package }


