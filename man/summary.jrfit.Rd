\name{summary.jrfit}
\alias{summary.jrfit}
\title{
Summary function for jrfit.
}
\description{
Returns a summary table for with estimates, standard errors, t-ratios, and p-values.
}
\usage{
\method{summary}{jrfit}(object,int=!object$fitint,df=NULL,...)
}
\arguments{
  \item{object}{The result of a call to jrfit }
  \item{int}{Should first row (intercept) of summary be printed? Defaults to FALSE if itercept was fit in \code{jrfit}.}
  \item{df}{Degrees of freedom to use in calculation of p-values.}
  \item{\dots}{additional arguments. currently unused. }

}

\details{
Forms standard table of estimates, standard errors, Wald type t-ratios and p-values.

Uses \code{df} = N - ncol(X) - 1 if \code{var.type} was set to 'cs' in \code{jrfit} where ncol(X) is the number of columns in the full design matrix including the columns for the block effects if there were fit too.  
Otherwise \code{df} = m.
}

\references{
Kloke, J.D., McKean, J.W., Rashid, M. (2009), Rank-based estimation and associated inferences for linear models with cluster correlated errors, \emph{Journal of the American Statistical Association}, 104, 384-390.
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\examples{
data(crabgrass)
x<-crabgrass[,1:6]
y<-crabgrass[,7]
block<-crabgrass[,8]
summary(jrfit(x,y,block))
}
