## wccStick.R

## access to the C functions in wccStick.c. Both pat1 and pat2 should
## be two-column matrices, with rt in the first column and I in the
## second. Returns the wcc value for these two stick spectra.
wcc.st <- function(pat1, pat2, trwidth) {
  WCC <- 0
  np1 <- nrow(pat1)
  np2 <- nrow(pat2)
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

## #######################################################################
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
  ## first we take out the mass info that is not used here
  ref <- lapply(ref, function(x) x[, c("rt", "I"), drop = FALSE])
  samp <- lapply(samp, function(x) x[, c("rt", "I"), drop = FALSE])
                 
  a <- init.coef
  if (!try) {
    Opt <- optim(a, STWCC, NULL, ref, samp, trwdth = trwdth, ...)
    
    a <- c(Opt$par)
  }

  v <- STWCC(a, ref, samp, trwdth.res)

  list(a = a, v = v)
}


STWCC <- function(warp.coef, refList, sampList, trwdth) {
  for (i in 1:length(sampList)) {
    sampList[[i]][,"rt"] <- warp.time(sampList[[i]][,"rt"], warp.coef)
    sampList[[i]] <-
        sampList[[i]][!is.na(sampList[[i]][,"rt"]), , drop = FALSE]
  }
  
  wccs <- mapply(wcc.st, refList, sampList, trwidth = trwdth)
  
  1 - mean(wccs)
}

warp.time <- function(tp, coef) {
  powers <- 1:length(coef) - 1

  ## the next correction is the wrong way around.
  new.tp <- c(outer(tp, powers, FUN = "^") %*% coef)
  ## We assume some form of symmetry and add the diffs with the
  ## original points to the list

  tp2 <- outer(2*tp - new.tp, -3:3, "*")
  tp2 <- sort(c(tp2))

  new.tp <- c(outer(tp2, powers, FUN = "^") %*% coef)
  approx(new.tp, tp2, xout = tp)$y
}
## warp.time <- function(tp, coef) {
##   powers <- 1:length(coef) - 1
##   ## it may happen that for an mz channel there is only one rt and
##   ## therefore the difference is zero: in that case we simply take
##   ## 5 seconds. If there is a range of rt values, we take 1 percent of
##   ## that range.  
##   small.value <- max(diff(tp[c(1, ntp)])*0.01, 1)

##   ## the next correction is the wrong way around.
##   new.tp <- c(outer(tp, powers, FUN = "^") %*% coef)
##   ## We assume some form of symmetry and add the diffs with the
##   ## origina points to the list
  
##   tp <- sort(tp)
##   ntp <- length(tp)

##   tp2 <- sort(c(2*tp[1] - tp[2],
##                 rbind(tp - small.value, tp + small.value),
##                 2*tp[ntp-1] - tp[ntp]))

##   new.tp <- c(outer(tp2, powers, FUN = "^") %*% coef)
##   approx(new.tp, tp2, xout = tp)$y
## }

## stick version of ptw: always global alignment, always using WCC, no
## selected traces,  no try argument. Here ref and sample are derived
## from peak tables as generated, e.g., by xcms, separated out into
## different mass channels - each channel is a list element. Function
## pktab2mzchannel is doing this.
stptw <- function (ref, samp, 
                   init.coef = c(0, 1, 0), 
                   trwdth = 20, trwdth.res = trwdth, ... )
{
  WCC <- stwarp(ref, samp, c(init.coef),
                trwdth = trwdth, trwdth.res = trwdth.res,
                ...)

  warped.sample <- lapply(samp,
                          function(x) {
                            x[,"rt"] <- warp.time(x[, "rt"], WCC$a)
                            x
                          })
  
  result <- list(reference = ref, sample = samp,
                 warped.sample = warped.sample,
                 warp.coef = WCC$a, warp.fun = NULL,
                 crit.value = WCC$v, optim.crit = "WCC",
                 warp.type = "global")
  class(result) <- c("stptw", "ptw")
  
  result
}


## ###################################################################
## reasonably brainless pieces of code

summary.stptw <- function (object, ...) {
  nsamp <- length(object$sample)
  nref <- length(object$reference)
  cat("PTW object:", object$warp.type,
      ifelse((object$warp.type == "individual" & nsamp > 1),
             "alignments of", "alignment of"), 
      nsamp, ifelse(nsamp > 1, "samples on", "sample on"), 
      nref, ifelse(nref > 1, "references.\n", "reference.\n"))
  cat("\nWarping coefficients:\n")
  print(object$warp.coef)
  cat("\nWarping criterion:", object$crit.type)
  cat(ifelse(object$warp.type == "individual" & nsamp > 1, 
             "\nValues:", "\nValue:"), object$crit.value, "\n\n")
}

print.stptw <- function (object, ...) {
  nsamp <- length(object$sample)
  nref <- length(object$reference)
  cat("PTW object:", object$warp.type,
      ifelse((object$warp.type == "individual" & nsamp > 1),
             "alignments of", "alignment of"), 
      nsamp, ifelse(nsamp > 1, "samples on", "sample on"), 
      nref, ifelse(nref > 1, "references.\n", "reference.\n"))
}

## ###################################################################
## Aux functions for converting peak tables into lists of mz channels

pktab2mzchannel <- function(pktab, Ivalue = "maxo", masses = NULL,
                            nMasses = 0, massDigits = 2) {
  pkst <- pktab[, c("mz", "rt", Ivalue)]
  colnames(pkst)[3] <- "I"

  if (is.null(masses)) {
    pkst[,"mz"] <- round(pkst[,"mz"], massDigits)

    if (nMasses == 0) {
      ## we assume they are already sorted!
      masses <- unique(pkst[,"mz"])
    } else {
      massTab <- sort(table(pkst[,"mz"]), decrease = TRUE)
      masses <- as.numeric(names(massTab))[1:nMasses]
      masses <- sort(masses)
    }    
  }

  threshold <- 5*10^{-massDigits-1}
  mzdiff <- abs(outer(pkst[,"mz"], masses, "-"))
  result <- lapply(1:length(masses),
                   function(ii)
                   pkst[which(mzdiff[,ii] < threshold),,drop=FALSE])
  names(result) <- masses

  result
}

mzchannel2pktab <- function(mzchannels) {
  do.call("rbind", mzchannels)
}

## ###################################################################

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
    

