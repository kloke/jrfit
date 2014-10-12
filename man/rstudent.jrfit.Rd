\name{rstudent.jrfit}
\alias{rstudent.jrfit}
\title{
Studentized residuals for the JR fit.
}
\description{
Computes studentized residuals based on \code{jrfit}.
Assumes compound symmetric errors.
}
\usage{
\method{rstudent}{jrfit}(model,...)
}
\arguments{
  \item{model}{object of class jrfit, the result of a call to \code{jrfit} }
  \item{\dots}{ additional arguments.  currently not used.  }
}
\details{
Estimates of the variance components are based on MAD of the prediction.  See Kloke, McKean, Rashid (2009) for details.  
}
\value{
A vector of length n containing the studentized residuals.
}
\references{
Kloke, J.D., McKean, J.W., Rashid, M. (2009), Rank-based estimation and associated inferences for linear models with cluster correlated errors, \emph{Journal of the American Statistical Association}, 104, 384-390.
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\seealso{
\code{\link{jrfit}}
}
\examples{
# a simple simulated example w/ 4 equal sized blocks #
x<-matrix(rnorm(40),ncol=2); y<-rnorm(20) ; g<-rep(1:4,each=5)
# rstudent(jrfit(x,y,g))
}
