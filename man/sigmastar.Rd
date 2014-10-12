\name{sigmastar}
\alias{sigmastar}
\title{Estimate of sigma(0).}
\description{
An estimate of sigma(0) needed for an estimate of the standard error of the estimate of intercept.
}
\usage{
sigmastar(ehat, block, p)
}
\arguments{ 
\item{ehat}{ n by 1 vector of residuals } 
\item{block}{ n by 1 vector denoting block membership }
\item{p}{ number of regression parameters (rank of design matrix)  }
}
\references{
Kloke, J.D., McKean, J.W., Rashid, M. (2009), Rank-based estimation and associated inferences for linear models with cluster correlated errors, \emph{Journal of the American Statistical Association}, 104, 384-390.
}
\author{
John Kloke <kloke@biostat.wisc.edu>
}

\examples{
## The function is currently defined as
function (ehat, block, p) 
{
    ublock <- unique(block)
    m <- length(ublock)
    nvec <- vector(m, mode = "numeric")
    signe <- sign(ehat)
    total <- 0
    for (k in 1:m) {
        signek <- signe[block == ublock[k]]
        nvec[k] <- length(signek)
        for (i in 1:(nvec[k] - 1)) {
            for (j in (i + 1):nvec[k]) {
                total <- total + signek[i] * signek[j]
            }
        }
    }
    M <- sum(choose(nvec, 2)) - p
    rhos <- total/M
    nstar <- sum(nvec * (nvec - 1))/sum(nvec)
    1 + nstar * rhos
  }
}
