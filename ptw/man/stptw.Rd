\name{stptw}
\alias{stptw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stptw(ref, samp, init.coef = c(0, 1, 0), trwdth = 20, trwdth.res = trwdth, nGlobal = ifelse(length(init.coef) > 3, 5, 0), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
%%     ~~Describe \code{ref} here~~
}
  \item{samp}{
%%     ~~Describe \code{samp} here~~
}
  \item{init.coef}{
%%     ~~Describe \code{init.coef} here~~
}
  \item{trwdth}{
%%     ~~Describe \code{trwdth} here~~
}
  \item{trwdth.res}{
%%     ~~Describe \code{trwdth.res} here~~
}
  \item{nGlobal}{
%%     ~~Describe \code{nGlobal} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ref, samp, init.coef = c(0, 1, 0), trwdth = 20, trwdth.res = trwdth, 
    nGlobal = ifelse(length(init.coef) > 3, 5, 0), ...) 
{
    WCC <- stwarp(ref, samp, init.coef, trwdth = trwdth, trwdth.res = trwdth.res, 
        nGlobal = nGlobal, ...)
    warped.sample <- lapply(samp, function(x) {
        if (nrow(x) > 0) 
            x[, "rt"] <- warp.time(x[, "rt"], WCC$a)
        x
    })
    result <- list(reference = ref, sample = samp, warped.sample = warped.sample, 
        warp.coef = WCC$a, warp.fun = NULL, crit.value = WCC$v, 
        optim.crit = "WCC", warp.type = "global")
    class(result) <- c("stptw", "ptw")
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
