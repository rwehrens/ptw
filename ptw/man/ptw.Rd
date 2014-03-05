\name{ptw}
\alias{ptw}
\alias{pmwarp}
\alias{print.ptw}
\alias{summary.ptw}
\title{Parametric Time Warping}
\description{The main function of the \code{ptw} package: it is a
  wrapper for the function \code{pmwarp}, which performs parametric time
  warping of one or more samples. Features in the samples are optimally
  aligned with features in the reference(s). One may align a single sample 
  to a single reference, several samples to a single reference, and several 
  samples to several references. In the latter case, the number of 
  references and samples should be equal. One may require that all samples 
  are warped with the same warping function, or one may allow individual
  warpings for all samples.}
\usage{
ptw(ref, samp, selected.traces,
    init.coef = c(0, 1, 0), try = FALSE,
    warp.type = c("individual", "global"),
    optim.crit = c("WCC", "RMS"),
    smooth.param = ifelse(try, 0, 1e05),
    trwdth = 20, trwdth.res = trwdth,
    verbose = FALSE, ...)
pmwarp(ref, samp, optim.crit, init.coef, try = FALSE,
       trwdth, trwdth.res, smooth.param, ...)
\method{summary}{ptw}(object, \dots)
\method{print}{ptw}(x, \dots)
}
\arguments{
  \item{ref}{reference. Either a vector (containing one reference
    signal) or a matrix (one reference per row). If
    more than one reference is specified, the number of reference
    signals must equal the number of sample signals.} 
  \item{samp}{sample. A vector (containing one sample signal) or a
    matrix (one sample per row).} 
  \item{selected.traces}{optional vector containing the row numbers to
    use from \code{ref} (if more than one reference signal is specified)
    and \code{samp}.}   
  \item{init.coef}{starting coefficients. The first number is the
    zeroth-order coefficient (i.e., a constant shift); further numbers
    indicate linear, quadratic, ... stretches. The default is to start
    from the identity warping using a quadratic function
    (\code{c(0, 1, 0)})}  
  \item{try}{if \code{try = TRUE}, \code{ptw} does not optimize the
    warping but returns a \code{ptw} object containing the warping for
    \code{init.coef}. Default: \code{FALSE}}
  \item{warp.type}{default is to treat samples and references as single
    entities and align them individually and independently. Using the 
    argument \code{warp.type = "global"} leads to one alignment 
    function; the samples are warped simultaneously to the 
    reference(s). Also see details}
  \item{optim.crit}{either \code{"WCC"} or \code{"RMS"}. In both cases,
    the optimal value of the alignment leads to a value of 0. For "WCC",
    this means that \code{1 - WCC} is optimized rather than \code{WCC}
    (where the optimal value equals 1)}
  \item{smooth.param}{smoothing parameter for smoothing the reference
    and sample when \code{optim.crit} equals \code{"RMS"}. If no
    smoothing is required, set this to 0. The default is to use
    smoothing in the optimization mode, and no smoothing otherwise}
  \item{trwdth}{the width of the triangle in the WCC criterion during
    the optimization, given as a number of data points. Default: 20}   
  \item{trwdth.res}{the width of the triangle in the WCC calculation in
    the calculation of the quality of the final result. Default: equal
    to \code{trwdth}}
  \item{verbose}{logical, default is \code{FALSE}. Whether to give
    output during the optimisation, which may be useful for large data
    sets}
  \item{...}{further arguments to optim}
  \item{x, object}{an object of class "ptw"}
}
\details{In the optimization mode (\code{try = FALSE}), the function
  optimizes the warping coefficients using the chosen criterion (either
  "WCC" or "RMS"). For "RMS", the data are smoothed before the
  optimization, but the quality of the final warping is measured on the
  unsmoothed data. For "WCC", the warping is performed using 
  \code{trwdth} as the triangle width, but the quality of the final 
  solution is measured using \code{trwdth.res}.
  
  If \code{try = TRUE} is used as an argument, the function does not
  start an optimization, but just calculates the warping for the given 
  warp function (\code{init.coef}); if \code{smooth.param} is larger 
  than zero for the RMS criterion, the RMS of the smoothed patterns is 
  calculated. The WCC criterion uses \code{trwidth.res} as the triangle 
  width in this case. 
  
  Five situations can be distinguished:
  \enumerate{
  \item One sample and one reference: this obviously leads to one warping
    function regardless of the setting of \code{warp.type}.
  \item Several samples, all warped to the same single reference, each with 
    its own warping function: this is the default behaviour
    (\code{warp.type = "individual"})
  \item Several samples, warped to an equal number of references
    (pair-wise), with their own warping functions: this is the default
    behaviour (\code{warp.type = "individual"})
  \item Several samples, warped to one reference, with one warping function 
    (\code{warp.type = "global"})
  \item Several samples, warped to an equal number of references
    (pair-wise), with one warping function 
    (\code{warp.type = "global"})
  }
}
\value{A list of class "ptw" containing:
  \item{reference}{the reference(s) used as input}
  \item{sample}{the sample(s) used as input}
  \item{warped.sample}{the warped sample}
  \item{warp.coef}{the warping coefficients}
  \item{warp.fun}{the warped indices}
  \item{crit.value}{the value of the chosen criterion, either "WCC" or "RMS"}
  \item{optim.crit}{the chosen criterion, either "WCC" or "RMS"}
  \item{warp.type}{the chosen type of warping, either "individual" or
    "global"}
}
\seealso{
  \code{\link{WCC}}, \code{\link{RMS}}, \code{\link{select.traces}}
}
\references{
  Eilers, P.H.C. (2004) "Parametric Time Warping", Analytical Chemistry, \bold{76} (2), 404 -- 411.
  
  Bloemberg, T.G., et al. (2010) "Improved parametric time warping for Proteomics", Chemometrics and Intelligent Laboratory Systems, \bold{104} (1), 65 -- 74.
}
\examples{
data(gaschrom)
ref <- gaschrom[1,]
samp <- gaschrom[16,]
gaschrom.ptw <- ptw(ref, samp)
summary(gaschrom.ptw)

gaschrom.ptw <- ptw(ref, samp, init.coef = c(0, 1, 0, 0))
summary(gaschrom.ptw)

\dontrun{
ref <- gaschrom[1,]
samp <- gaschrom[2:16,]
gaschrom.ptw <- ptw(ref, samp, warp.type = "individual", verbose = TRUE,
              optim.crit = "RMS", init.coef = c(0, 1, 0, 0))
summary(gaschrom.ptw)

ref <- gaschrom[1:8,]
samp <- gaschrom[9:16,]
gaschrom.ptw <- ptw(ref, samp, warp.type = "individual",
              optim.crit = "RMS", init.coef = c(0, 1, 0, 0))
summary(gaschrom.ptw)

gaschrom.ptw <- ptw(ref, samp, warp.type = "global",
              optim.crit = "RMS", init.coef = c(0, 1, 0, 0))
summary(gaschrom.ptw)

# Example of a three-way data set
data(lcms)
# first bring all samples to the same scale
lcms.scaled <- aperm(apply(lcms, c(1,3), 
                           function(x) x/mean(x) ), c(2,1,3))
# add zeros to the start and end of the chromatograms
lcms.s.z <- aperm(apply(lcms.scaled, c(1,3), 
                        function(x) padzeros(x, 250) ), c(2,1,3))
# define a global 2nd degree warping
warp1 <- ptw(lcms.s.z[,,2], lcms.s.z[,,3], warp.type="global")
warp.samp <- warp1$warped.sample
warp.samp[is.na(warp.samp)] <- 0
# refine by adding 5th degree warpings for individual chromatograms
warp2 <- ptw(lcms.s.z[,,2], warp.samp, init.coef=c(0,1,0,0,0,0))
warp.samp2 <- warp2$warped.sample
warp.samp2[is.na(warp.samp2)] <- 0
# compare TICs
layout(matrix(1:2,2,1, byrow=TRUE))
plot(colSums(lcms.s.z[,,2]), type="l", ylab = "",
     main = "TIC: original data")
lines(colSums(lcms.s.z[,,3]), col=2, lty=2)
plot(colSums(lcms.s.z[,,2]), type="l", ylab = "",
     main = "TIC: warped data")
lines(colSums(warp.samp2), lty=2, col=2)
}
}
\author{Jan Gerretzen, Ron Wehrens}
\keyword{manip}
