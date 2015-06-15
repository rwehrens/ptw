\name{wcc.st}
\alias{wcc.st}
\alias{wac.st}
\title{
  Weighted auto- and cross-correlation measures for peak lists
}
\description{Functions to calculate weighted auto- and cross-correlation
  measures for peak lists; these functions are the analogues of
  functions \code{wcc} and \code{wac} for profiles. The WCC is a
  suitable measure for the similarity of two patterns when features may
  be shifted. Identical patterns lead to a WCC value of 1. 
}
\usage{
wcc.st(pat1, pat2, trwidth, acors1 = NULL, acors2 = NULL)
wac.st(pat1, trwidth)
}
\arguments{
  \item{pat1, pat2}{input patterns. Matrices of two columns, with
    retention time in the first column and intensity in the second.}
  \item{trwidth}{triangle width, given in the unit of time.}
  \item{acors1, acors2}{autocorrelations of the input patterns. If not
    provided, they are calculated}
}
\value{One number, the weighted autocorrelation or crosscorrelation}
\references{
  de Gelder, R., Wehrens, R. and Hageman, J.A. (2001) "A generalized expression for the similarity of spectra: Application to powder diffraction pattern classification",  Journal of Computational Chemistry, \bold{22}, 273 -- 289.

Wehrens, R., T.G. Bloemberg, T.G., and Eilers, P.H.C. (2015)
"Fast parametric warping of peak lists", Bioinformatics, in press. DOI: 10.1093/bioinformatics/btv299}
\author{Ron Wehrens}
\examples{
}
\keyword{manip}
