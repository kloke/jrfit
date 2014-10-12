\name{wald.test.jrfit}
\alias{wald.test.jrfit}
\title{
Do a Wald test based on a JR fit.
}
\description{
Do a Wald test based on a JR fit.
}
\usage{
wald.test.jrfit(fit, K, asymptotic = FALSE)
}
\arguments{
  \item{fit}{
Output from jrfit.
}
  \item{K}{
A full row rank contrast matrix.
}
  \item{asymptotic}{
Should asymptotic (chi-sq) values be used or small sample (F).
}
}
\details{
Test of K beta = 0.
}
\value{
\item{statistic}{The test statistic}
\item{p.value}{The p-value}
\item{asymptotic}{logical indicating if asymptotic critical values were used}
\item{df}{Which chi-sq or F degress of freedom were used}
}
\references{
Kloke, J.D., McKean, J.W., Rashid, M. (2009), Rank-based estimation and associated inferences for linear models with cluster correlated errors, \emph{Journal of the American Statistical Association}, 104, 384-390.
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\seealso{
\code{\link{jrfit}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fit, K, asymptotic = FALSE) 
{
    est <- K \%*\% fit$coef
    var.est <- K \%*\% fit$varhat \%*\% t(K)
    statistic <- t(est) \%*\% solve(var.est) \%*\% est
    q <- nrow(K)
    if (asymptotic) {
        p.value <- pchisq(statistic, q, lower.tail = FALSE)
        df <- q
    }
    else {
        statistic <- statistic/q
        p.value <- pf(statistic, q, fit$DF, lower.tail = FALSE)
        df <- c(q, fit$DF)
    }
    list(statistic = statistic, p.value = p.value, asymptotic = asymptotic, 
        df = df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
