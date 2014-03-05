\name{interpol}
\alias{interpol}
\title{Linear interpolation}
\description{This function calculates values at non-integer index
  positions by linear interpolation.}
\usage{interpol(xout, y)}
\arguments{
  \item{xout}{points at which to compute interpolated values}
  \item{y}{signal to be interpolated}
}
\value{A vector of interpolated values. If any of the values of
  \code{xout} are outside the domain \code{1:length(y)} the output will
  contain NA at those positions
}
\details{The C function approx.c that is called by interpol, was copied 
  from the \code{stats} package of R-2.15.1. This was done to prevent aberrant 
  behaviour due to changes without notice.
}
\references{R Development Core Team (2012). R: A language and environment for
  statistical computing. R Foundation for Statistical Computing,
  Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org/.
} 
\author{R Development Core Team (see details)}
\examples{
interpol(seq(23, 25, length = 20),10:50)
interpol(seq(23, 25, length = 20),5:50)
}
\keyword{manip}
