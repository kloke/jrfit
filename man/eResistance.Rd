\name{eResistance}
\alias{eResistance}
\docType{data}
\title{
Electrical Resistance Data
}
\description{
Results of an experiment designed to determine if five electrode types performed similarly.
Each electrode type (etype) was applied to the arm of 16 subjects.
The outcome variable is resistance.
}
\usage{data(eResistance)}
\format{
A matrix with 80 observations on three variables.
}
\source{
Stokes, Davis, Koch (1995), \emph{Categorical Data Analysis Using the SAS System}, Cary, NC: SAS Institute Inc.
}
\examples{
data(eResistance)
boxplot(resistance~etype,data=eResistance)
}
\keyword{datasets}
