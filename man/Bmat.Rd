\name{Bmat}
\alias{Bmat}
\title{ An internal function for jrfit}
\description{ An internal function for jrfit.  }
\usage{
Bmat(v1, v2, block)
}
\arguments{
  \item{v1}{scalar}
  \item{v2}{scalar}
  \item{block}{vector}
}
\author{ John Kloke <kloke@biostat.wisc.edu> }
\examples{
function (v1, v2, block) 
{
    ublock <- unique(block)
    m <- length(ublock)
    N <- length(block)
    B <- matrix(0, nrow = N, ncol = N)
    ind <- 1
    for (k in 1:m) {
        nk <- sum(block == ublock[k])
        Bk <- diag(rep(v1 - v2, nk)) + v2
        B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- Bk
        ind <- ind + nk
    }
    B
  }
}
