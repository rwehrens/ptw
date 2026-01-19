\docType{package}
\name{ptw-package}
\alias{ptw-package}
\title{ptw: Parametric Time Warping}
\description{
Parametric Time Warping aligns patterns, i.e. it aims to put corresponding features at the same locations. The algorithm searches for an optimal polynomial describing the warping. It is possible to align one sample to a reference, several samples to the same reference, or several samples to several references. One can choose between calculating individual warpings, or one global warping for a set of samples and one reference. Two optimization criteria are implemented: RMS (Root Mean Square error) and WCC (Weighted Cross Correlation). Both warping of peak profiles and of peak lists are supported. A vignette for the latter is contained in the inst/doc directory of the source package - the vignette source can be found on the package github site. See `citation("ptw")` for more details.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/rwehrens/ptw}
  \item Report bugs at \url{https://github.com/rwehrens/ptw/issues}
}

}
\author{
\strong{Maintainer}: Ron Wehrens \email{ron.wehrens@gmail.com}

Authors:
\itemize{
  \item Paul Eilers
  \item Tom Bloemberg \email{t.bloemberg@science.ru.nl}
}

Other contributors:
\itemize{
  \item Jan Gerretzen [contributor]
  \item Hans Wouters [contributor]
}

}
\keyword{internal}
