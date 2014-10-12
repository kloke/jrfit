\name{print.jrfit}
\alias{print.jrfit}
\alias{print.summary.jrfit}
\title{jrfit Internal Print Functions}
\description{
These functions print the output in a user-friendly manner using the internal R function \code{print}.
}
\usage{
\method{print}{jrfit}(x, ...)
\method{print}{summary.jrfit}(x, digits = max(5, .Options$digits - 2), ...)
}
\arguments{
  \item{x}{ An object to be printed }
  \item{digits}{number of digits to display}
  \item{\dots}{  additional arguments to be passed to \code{print} }
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\seealso{
\code{\link{jrfit}}, \code{\link{summary.jrfit}}
}
