\name{jrfit-package}
\alias{jrfit-package}
%\alias{jrfit}
\docType{package}
\title{Joint/Jaeckel Rankings (JR) estimates and inference for cluster correlated data.}
\description{
JR estimates were developed by Kloke, et. al. (2009) and extend the results Jureckova (1971) and Jaeckel (1972) to estimating parameters in a mixed model with one random effect.
Inference in complete including tests of estimates, standard errors, tests of hypothesis, studentized residuals.
}
\details{
\tabular{ll}{
Package: \tab jrfit\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-10-10\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
}
\author{
John Kloke

Maintainer: John Kloke <kloke@biostat.wisc.edu>
}
\references{
Jaeckel, L.A. (1972), Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annal
s of Mathematical Statistics}, 43, 1449 - 1458.

Jureckova, J. (1971), Nonparametric estimate of regression coefficients. \emph{Annals of Mathematical Statistics
}, 42, 1328 - 1338.

Kloke, J.D., McKean, J.W., Rashid, M. (2009), Rank-based estimation and associated inferences for linear models with cluster correlated errors, \emph{Journal of the American Statistical Association}, 104, 384-390.

}
%\keyword{ package }
%\seealso{
%\code{\link[Rfit:<pkg>-package]{<pkg>}} ~~
%}
\examples{
data(crabgrass)
x<-crabgrass[,1:6]
y<-crabgrass[,7]
block<-crabgrass[,8]
jrfit(x,y,block)
}
