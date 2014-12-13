\name{vee}
\alias{vee}
\title{Estimates of the Variance Components}
\description{
For cluster correlated data, 
estimates of the variance components are computed as discussed in Section 8.4 of Kloke and McKean (2014).
There are two choices: if the value of method is mm then the medain and MAD estimators while
if it has the value dhl the Hodges-Lehmann estimate and the Wilcoxon dispersion function is used.
The mm type are robust estimators while the dhl type are efficient estimators.
}
\usage{
vee(ehat,center,method='dhl',scores=wscores)
}
\arguments{
  \item{ehat}{ N x 1 vector of residuals }
  \item{center}{ N x 1 vector denoteing the center or cluster }
  \item{method}{ option indicator as discussed above }
  \item{scores}{ score function used in the fit of ehat }
}
\details{
Estimates are discussed in detail in Section 8.4 of Kloke and McKean (2014).
}
\value{
   \item{sigb2}{estimate of the variance of the random effect}
   \item{sige2}{estimate of the variance of the random error}
}
\references{
Kloke, J.D., McKean, J.W., and Rashid, M. (2009), Rank-based estimation and associat
ed inferences for linear models with cluster correlated errors, \emph{Journal of
 the American Statistical Association}, 104, 384-390.

Kloke, J. and McKean, J.W. (2014), \emph{Nonparametric statistical methods using R}, Boca Raton, FL: Chapman-Hall.

}
\author{ John Kloke <kloke@biostat.wisc.edu> }

\examples{

m<-10  # number of blocks
n<-10  # number number
k<-2   # number of treatments
N<-m*n # total sample size

x<-rnorm(N)                       # covariate
w<-sample(c(0,1),N, replace=TRUE) # treatment indicator
block<-rep(1:m,n)                 # m blocks of size n

X<-cbind(x,w)

Z<-model.matrix(~as.factor(block)-1)
b<-rnorm(m,sd=3)
e<-rnorm(N)
y<-Z\%*\% b+e

fit<-jrfit(X,y,block)
summary(fit)
vee(fit$resid,fit$block,method='mm')

vee(fit$resid,fit$block)
}

