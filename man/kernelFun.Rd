\name{kernelFun}
\alias{kernelFun}

\title{
Several kernel functions
}
\description{
Computes several kernel functions(truncated, Bartlett, Daniell, QS, Parzen).
These kernels are for constructing test statistics for testing pairwise
independence.
}
\usage{
kernelFun(type, z)
}

\arguments{
\item{type}{
A character string which indicates the name of the smoothing kernel.
\code{kernelFun} can be: 'truncated', 'bartlett', 'daniell', 'QS', 'parzen'.
No default is given.
}
\item{z}{
A real number.
}
}
\details{
\code{kernelFun} computes several kernel functions including truncated,
Bartlett, Daniell, QS and Parzen.

The exact definition of each of the above functions are given below: \cr
\itemize{
\item Truncated
\deqn{ k(z) = \left\{
                     \begin{array}{ll}
                       1, & |z| \leq 1, \\[1ex]
                       0, & \mbox{otherwise}.
                     \end{array}
                   \right. }
\item Bartlett
\deqn{ k(z) = \left\{
                     \begin{array}{ll}
                       1 - |z|, & |z| \leq 1, \\[1ex]
                       0, & \mbox{otherwise}.
                     \end{array}
                   \right. }
\item Daniell
\deqn{ k(z) = \frac{\mbox{sin}(\pi z)}{\pi z}, z \in \Re - \{0\} }
\item QS
\deqn{ k(z)=(9/5\pi^2z^2)\{\mbox{sin}(\sqrt{5/3}\pi z)/\sqrt{5/3}\pi z-\mbox{cos}(\sqrt{5/3}\pi z)\}, z \in \Re }
\item Parzen
\deqn{ k(z) = \left\{
                          \begin{array}{ll}
                            1-6(\pi z/6)^2 + 6|\pi z/6|^3, & |z| \leq 3/\pi, \\[1ex]
                            2(1-|\pi z/6|)^3, & 3/\pi \leq |z| \leq 6/\pi, \\[1ex]
                            0, & \mbox{otherwise}
                          \end{array}
                        \right. }
}
All these kernel functions are mainly used to smooth the generalized spectral
density function, firstly introduced by Hong (1999). Assumptions and theoretical
properties of these functions can be found in Hong (1996;1999) and
Fokianos and Pitsillou (2017).
}

\value{
A value that lies in the interval \eqn{[-1, 1]}.
}
\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature
Review of Distance Correlation and Its Applications to Time Series.
\emph{International Statistical Review}, 87, 237-262.

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise dependence
in time series. \emph{Technometrics}, 159, 262-3270.

Pitsillou M. and Fokianos K. (2016). dCovTS: Distance Covariance/Correlation
for Time Series. \emph{R Journal}, 8, 324-340.

Hong, Y. (1996). Consistent testing for serial correlation of unknown form.
\emph{Econometrica}, 64, 837-864.

Hong, Y. (1999). Hypothesis testing in time series via the empirical
characteristic function: A generalized spectral density approach.
\emph{Journal of the American Statistical Association}, 94, 1201-1220.
}
\author{
Maria Pitsillou and Konstantinos Fokianos.
}
\examples{
k1 <- kernelFun( "bartlett", z = 1/3 )
k2 <- kernelFun( "bar", z = 1/5 )
k3 <- kernelFun( "dan", z = 0.5 )
}

\keyword{ ts }

