\name{kernelFun}
\alias{kernelFun}

\title{
Several kernel functions
}
\description{
Computes several kernel functions(truncated, Bartlett, Daniell, QS, Parzen). These kernels are
for constructing test statistics for testing pairwise independence.
}
\usage{
kernelFun(type, z)
}

\arguments{
  \item{type}{
character string which indicates the name of the smoothing kernel. 
\code{kernelFun} can be: 'truncated', 'bartlett', 'daniell', 'QS', 'parzen'. No default is given.
}
  \item{z}{
real number.
}
}
\details{
\code{kernelFun} computes several kernel functions including truncated, Bartlett, Daniell, QS and Parzen.

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
All these kernel functions are mainly used to smooth the generalized spectral density function, firstly introduced by Hong (1999).
Assumptions and theoretical properties of these functions can be found in Hong (1996;1999) and Fokianos and Pitsillou (2016).
}

\value{
A value that lies in the interval \eqn{[-1,1]}.
}
\references{
Fokianos K. and M. Pitsillou (2016). Consistent testing for pairwise dependence in time series. \emph{Technometrics},
\url{http://dx.doi.org/10.1080/00401706.2016.1156024}.

Hong, Y. (1996). Consistent testing for serial correlation of unknown form. \emph{Econometrica} \eqn{\textbf{64}}{64}, 837-864,
\url{http://dx.doi.org/10.2307/2171847}.

Hong, Y. (1999). Hypothesis testing in time series via the empirical characteristic function: A generalized
spectral density approach. \emph{Journal of the American Statistical Association} \eqn{\textbf{94}}{94}, 1201-1220,
\url{http://dx.doi.org/10.1080/01621459.1999.10473874}.
}
\author{
Maria Pitsillou and Konstantinos Fokianos
}
\examples{
k1 <- kernelFun("bartlett",z=1/3)
k2 <- kernelFun("bar",z=1/5)
k3 <- kernelFun("dan",z=0.5)
}

\keyword{ ts }

