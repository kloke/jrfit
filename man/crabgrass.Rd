\name{crabgrass}
\alias{crabgrass}
\docType{data}
\title{
Crabgrass data.
}
\description{
Data are presented in Cobb (1998).
Split plot design (with two plots).  
Outcome variable is yield of crab grass.  
}
\usage{data(crabgrass)}
\format{
A matrix frame with 64 observations on 8 variables.
}
\details{
Cobb (1998) presented an example of a complete block
design concerning the weight of crab grass. The fixed factors in
the experiment are the density of the crabgrass (four levels) and
the levels (two) of the three nutrients nitrogen, phosphorus, and
potassium. Two complete blocks of the experiment were carried out, so altogether there are N = 64 observations. 
The response is the dry weight crab grass was grown in a cup.
}
\source{
Cobb (1998), \emph{Introduction to Design and Analysis of Experiments}, New York: Wiley.
}
\references{
Kloke, J.D., McKean, J.W., Rashid, M. (2009), Rank-based estimation and associated inferences for linear models with cluster correlated errors, \emph{Journal of the American Statistical Association}, 104, 384-390.
}
\examples{
data(crabgrass)
x<-crabgrass[,1:6]
y<-crabgrass[,7]
block<-crabgrass[,8]
summary(jrfit(x,y,block))
}
\keyword{datasets}
