## wccStick.R

## full R implementation of stick-based WCC for validation purposes
## and speed comparisons.

## R implementation of Cfg
## pat1 and pat2 two-column matrices (rt, I). Variable rtwidth is in
## the same scale as the rt. Eventually this function should be called
## in a loop for many different mz bins.

Rwght <- function(dff, trwidth) {
  1 - dff/trwidth
}

RCfg <- function(pat1, pat2, trwidth) {
  rtdiffs <- abs(outer(pat1[,1], pat2[,1], "-"))
  intr <- which(rtdiffs < trwidth, arr.ind = TRUE)
  if (nrow(intr) == 0) return(0)

  wghts <- sapply(1:nrow(intr),
                  function(ii)
                  Rwght(rtdiffs[ intr[ii,1], intr[ii,2] ],
                        trwidth))
  prods <- sapply(1:nrow(intr),
                  function(ii)
                  pat1[intr[ii,1],2] * pat2[intr[ii,2],2])
  sum(wghts*prods)
}

Rwcc.st <- function(pat1, pat2, trwidth) {
  RCfg(pat1, pat2, trwidth) /
      sqrt(RCfg(pat1, pat1, trwidth) * RCfg(pat2, pat2, trwidth))
}
    

## access to the C functions in wccStick.c

wcc.st <- function(pat1, pat2, trwidth) {
  WCC <- 0
  np1 <- nrow(pat1)
  np2 <- nrow(pat2)
  ##   cat("\nNp1:", np1, "Np2:", np2, "\n")
  p1 <- c(pat1)
  p2 <- c(pat2)
  
  res <- .C("st_WCC",
            as.double(p1),
            as.integer(np1),
            as.double(p2),
            as.integer(np2),
            as.double(trwidth),
            as.double(WCC),
            package = "ptw")
  
  res[[6]]
}


## Function to calculate a warped time point, given warping
## coefficients. This looks so easy... but it corrects exactly the
## wrong way. Hmm. Temporary fix...
warp.time <- function(tp, coef) {
  powers <- 1:length(coef) - 1
  new.tp <- c(outer(tp, powers, FUN = "^") %*% coef)

  2*tp - new.tp
}

## analogon to pmwarp.R: the actual optimization. Changes w.r.t. pmwarp:
## 1) only present for WCC, so all references to optim.crit are taken out
## 2) call optim using the STWCC function rather than WCC
## 3) for this initial version, ref.acors is not used, is calculated in C
## 4) for this initial version, wghts is not used, is calculated in C
## 5) trwdth and trwdth.res should be given in REAL TIME UNITS, and
##    not in terms of sampling points
stwarp <- function (ref, samp, init.coef, try = FALSE, trwdth, 
                    trwdth.res = trwdth, ...) 
{
  a <- init.coef
  if (!try) {
    Opt <- optim(a, STWCC, NULL, ref, samp, trwdth = trwdth, ...)
    
    a <- c(Opt$par)
  }

  v <- STWCC(a, ref, samp, trwdth.res)

  list(a = a, v = v)
}


## STWCC is the stick equivalent to WCC for peak-picked data,
## specifically designed for LC-MS data. This function gives the WCC
## value for a comparison between the reference (unchanged) and a
## particular warping of the sample. Both refList and sampList are lists
## of -eventually- three-column (mz, rt, I) matrices. For the moment
## we leave out the mass information. We could be aligning TIC values,
## so the matrices have only two columns (rt, I). Eventually, this
## should be called for all individual mz values.

STWCC <- function(warp.coef, refList, sampList, trwdth) {
  ## find new time points
  for (i in 1:length(sampList))
      sampList[[i]][,"rt"] <- warp.time(sampList[[i]][,"rt"], warp.coef)
  
  wccs <- mapply(wcc.st, refList, sampList, trwidth = trwdth)
  
  1 - mean(wccs)
}
