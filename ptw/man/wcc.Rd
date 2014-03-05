\name{wcc}
\alias{wcc}
\alias{wac}
\title{
  Weighted auto- and cross-correlation measures
}
\description{Functions to calculate weighted auto- and cross-correlation
  measures. The wcc is a suitable measure for the similarity of two
  patterns when features may be shifted. Identical patterns lead to a
  wcc value of 1.
}
\usage{
wcc(pattern1, pattern2, trwdth, wghts = NULL, acors1 = NULL, acors2 = NULL)
wac (pattern1, trwdth, wghts = NULL)
}
\arguments{
  \item{pattern1, pattern2}{input patterns, typically spectra. Vectors}
  \item{trwdth}{triangle width, given in the number of data points}
  \item{wghts}{optional weights vector, will be calculated from triangle
    width if necessary. Sometimes it is more efficient to
    pre-calculate it and give it as an argument}
  \item{acors1, acors2}{autocorrelations of the input patterns. If not
    provided, they are calculated}
}
\details{wcc and wac are defined such that the triangle width stands for 
  the number of points on one side of \strong{and including} the current 
  point. Thus, a \code{trwdth} of 0 signifies a non-existent triangle and 
  results in an error; a \code{trwdth} equal to 1 only includes the current
  point with weight 1 and no neighbouring points.

}
\value{One number, the weighted autocorrelation or crosscorrelation}
\references{
  de Gelder, R., Wehrens, R. and Hageman, J.A. (2001) "A generalized expression for the similarity of spectra: Application to powder diffraction pattern classification",  Journal of Computational Chemistry, \bold{22}, 273 -- 289.
}
\author{Ron Wehrens}
\examples{
data(gaschrom)
wcc(gaschrom[1,], gaschrom[2,], 20)
}
\keyword{classif}
