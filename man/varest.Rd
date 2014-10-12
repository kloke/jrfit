\name{varest}
\alias{tcs}
\alias{tsand}
\alias{tE}
\alias{tind}
\title{Estimates of the covariance.}
\description{
Estimates of the variance covariance of the scores.
Used in estimating variance  covariance of the linear model estimates and also in the score test.
}
\usage{
tcs(ehat, X, block, ahat = NULL, scores = wscores, eps = 1e-04)
}
\arguments{
  \item{ehat}{ N x 1 vector of residuals }
  \item{X}{ N x p design matrix }
  \item{block}{ N x 1 vector denoting block membership }
  \item{ahat}{ optional vector of scores evaluated at the residual (if null then is computed using ehat and scores) }
  \item{scores}{ score function used in the fit of ehat }
  \item{eps}{ used in tcs to ensure positive definite }
}
\details{
Estimates vhat =  sum_k^m X_k^T varphi_k X_k.  See Kloke, McKean, Rashid (2009) for details.
}
\value{
p x p matrix 
}
\references{
Kloke, J.D., McKean, J.W., Rashid, M. (2009), Rank-based estimation and associat
ed inferences for linear models with cluster correlated errors, \emph{Journal of
 the American Statistical Association}, 104, 384-390.
}
\author{ John Kloke <kloke@biostat.wisc.edu> }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ehat, X, block, ahat = NULL, scores = wscores, eps = 1e-04) 
{
    if (is.null(ahat)) 
        ahat <- getScores(scores, ehat)
    nvec <- tapply(ehat, block, length)
    M <- sum(choose(nvec, 2)) - ncol(X)
    myfunc <- function(ak) {
        Ak <- tcrossprod(ak)
        Ak[upper.tri(Ak)]
    }
    rhohat <- sum(unlist(tapply(ahat, block, myfunc)))/M
    K0 <- -1/(max(nvec) - 1)
    if (rhohat < K0) 
        rhohat <- K0 + eps
    if (rhohat > 1) 
        rhohat <- 1 - eps
    myfunc <- function(i, x, rhohat) {
        xi <- as.matrix(x[i, ])
        ni <- nrow(xi)
        sk <- matrix(rhohat, nrow = ni, ncol = ni)
        diag(sk) <- 1
        t(xi) \%*\% sk \%*\% xi
    }
    Reduce("+", tapply(1:nrow(X), block, myfunc, X, rhohat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
