\name{stptw}
\alias{stptw}
\alias{print.stptw}
\alias{summary.stptw}
\title{Parametric warping of peak lists}
\description{Where \code{ptw} is used for warping complete time
  profiles, this function is meant to align lists of peaks, containing
  retention time and intensity information. Peaks may be defined in
  different channels: e.g., in mass spectrometry each channel may be a
  specific m/z range; and in UV-Vis spectroscopy a channel may be a
  wavelength range or a component from Multivariate Curve
  Resolution. Compared to the \code{ptw} function, \code{stptw} has
  fewer options: it always uses \'forward\' warping, 
  the warp.type is \'global\', and the distance function always is \'WCC\'.
}
\usage{
stptw(ref, samp, init.coef = c(0, 1, 0), trwdth = 20, trwdth.res = trwdth, nGlobal = ifelse(length(init.coef) > 3, 5, 0), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{reference. A list of matrices, where each matrix is a set
    of peaks. Different list elements may refer to peaks measured at
    different wavelengths or different m/z values. The matrices do not
    need to have the same number of rows. Columns need to have names
    \"rt\" and \"I\".
}
  \item{samp}{Sample. A list of matrices, where each matrix is a set
    of peaks. The length of the list should be the same as the length of
    \code{ref}. Different list elements may refer to peaks measured at
    different wavelengths or different m/z values. The matrices do not
    need to have the same number of rows. Columns need to have names
    \"rt\" and \"I\".}
  \item{init.coef}{starting coefficients. The first number is the
    zeroth-order coefficient (i.e., a constant shift); further numbers
    indicate linear, quadratic, ... stretches. The default is to start
    from the identity warping using a quadratic function
    (\code{c(0, 1, 0)})}
  \item{trwdth}{the width of the triangle in the WCC criterion during
    the optimization, given as a number of data points. Default: 20}   
  \item{trwdth.res}{the width of the triangle in the WCC calculation in
    the calculation of the quality of the final result. Default: equal
    to \code{trwdth}}
  \item{nGlobal}{number of restarts for global optimisation. By default
    this is set to zero for quadratic or simpler warpings; for cubic or
    higher warpings the default is 5.}
  \item{\dots}{additional arguments to the optim}
}
\details{Apart from the differences mentioned above, the complete
  functionality of the \code{ptw} function is duplicated in
  \code{stptw}.

  Global optimisation using tools from the \code{nloptim} package is
  used to increase the chance of finding a global optimum. We have
  observed that for higher-order warping coefficients one sometimes can
  find very fluctuating results - the price of the improvement in
  stability is a much longer calculation time. This feature is not
  available in the profile-based alignment since calculation times would
  be too long, but with the improved speed in peak-based alignment it
  becomes feasible.}
\value{A list of class "stptw" containing:
  \item{reference}{the reference(s) used as input}
  \item{sample}{the sample(s) used as input}
  \item{warped.sample}{the warped sample}
  \item{warp.coef}{the warping coefficients}
  \item{warp.fun}{the warped indices}
  \item{crit.value}{the value of the chosen criterion, either "WCC" or
  "RMS". Included for compatibility with the "ptw" class}
  \item{optim.crit}{the chosen criterion, either "WCC" or
  "RMS". Included for compatibility with the "ptw" class} 
  \item{warp.type}{the chosen type of warping, either "individual" or
    "global". Included for compatibility with the "ptw" class}
}
\references{
R. Wehrens, T.G. Bloemberg and P.H.C Eilers. Fast parametric warping of
  peak lists. Bioinformatics (2015), in press. DOI:
  10.1093/bioinformatics/btv299}
\author{Ron Wehrens}
\seealso{
  \code{\link{ptw}}
}
\examples{
}
\keyword{manip}
