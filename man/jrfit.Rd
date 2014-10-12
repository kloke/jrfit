\name{jrfit}
\alias{jrfit}
\title{ Joint ranking estimates for cluster correlated data.  }
\description{
Returns R estimates for the linear models with cluster correlated errors.
Also returns objects useful for inference and diagnostics.
}
\usage{ jrfit(x, y, block, yhat0 = NULL, scores = wscores, fitint = NULL, 
var.type = 'sandwich',fitblock=FALSE,tuser=NULL,...) }
\arguments{
  \item{x}{ N by p design matrix }
  \item{y}{ N by 1 response vector }
  \item{block}{ N by 1 vector denoting cluster (block) membership }
  \item{yhat0}{ initial fit (defaults to computing l1 fit) }
  \item{scores}{ which score function to use (defaults to wscore (Wilcoxon scores) }
  \item{fitint}{ should an intercept be fit in addtion to the regression parameters (adds a col of 1s to x) }
  \item{var.type}{ one of sandwich (default), cs for compound symmetry, or user for user defined}
  \item{fitblock}{ should blocks be fit as (nuisance) fixed effects}
  \item{tuser}{ optional function to compute V_varphi. used when var.type='user' }
  \item{\dots}{ additional arguments. currently unused. }

}
\details{
Solves the rank based minimization problem using Jaeckel's (1972) dispersion function.
That is the ranks are taken over the entire dataset.
Results are presented in Kloke, et. al. (2009).

\code{scores} are available in code{\link{rfit}}.

If \code{fitint} is set to TRUE, a column of ones is added to the design matrix.  If it is set to NULL then, if 1 is not in the column space of the \code{x} then a column of ones is added.

The default behavior for \code{fitblock} is to set to TRUE when \code{var.type} is set to 'cs' and FALSE otherwise.
It is not recommended to set \code{fitblock} to TRUE and use \code{var.type} as 'sandwich'.
Setting \code{fitblock} to FALSE when using \code{var.type} as 'cs' may be useful at times.

\code{var.type} specifies how the variance covariance matrix of the parameter estimates should be estimated.
The default is to use a sandwich estimate.  Another option we have developed is the compound symmetry estimate (see Kloke, et. al. 2009).  The user is welcome to supply his or her own variance covariance function using the option 'user'.
This requires the user define \code{tuser}.
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449 - 1458.

Jureckova, J. (1971). Nonparametric estimate of regression coefficients. \emph{Annals of Mathematical Statistics}, 42, 1328 - 1338.

Kloke, J.D., McKean, J.W., Rashid, M.M. (2009).  Rank-based estimation and associated inferences for linear models with cluster correlated errors.  \emph{Journal of the American Statistical Association}, 104, 384-390.
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\seealso{
\code{\link{rfit}}
}
\examples{
# a simple simulated example w/ 4 equal sized blocks #
x<-matrix(rnorm(40),ncol=2); y<-rnorm(20) ; g<-rep(1:4,each=5)
jrfit(x,y,g)
}
\keyword{robust}
\keyword{nonparametric }
