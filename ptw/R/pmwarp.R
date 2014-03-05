pmwarp <- function (ref, samp, optim.crit, init.coef, try = FALSE,
                    trwdth, trwdth.res, smooth.param, ...)
{
  ## Multiply coefficients to prevent them from becoming too 
  ## small for numeric precision.
  n <- length(init.coef)
  ncr <- ncol(ref)
  time <- (1:ncr) / ncr
  B <- matrix(time, nrow = ncr, ncol = n)
  B <- t(apply(B, 1, cumprod))/B
  a <- init.coef * ncr^(0:(n-1))

  if (!try) { # perform optimization
    switch(optim.crit,
           RMS = {
             if (smooth.param > 0) {
               samp.sm <- t(apply(samp, 1, difsm, smooth.param))
               ref.sm <- t(apply(ref, 1, difsm, smooth.param))
               Opt <- optim(a, RMS, NULL, ref.sm, samp.sm, B, ...)
             } else {
               Opt <- optim(a, RMS, NULL, ref, samp, B, ...)
             }},
           WCC = {
             wghts <- 1 - (0:trwdth)/trwdth
             ref.acors <- apply(ref, 1, wac, trwdth = trwdth, wghts = wghts)
             Opt <- optim(a, WCC, NULL, ref, samp, B,
                          trwdth = trwdth, wghts = wghts,
                          ref.acors = ref.acors, ...)
           })
    
    a <- c(Opt$par)
    v <- Opt$value

    if ((optim.crit == "RMS" && smooth.param > 0) ||
        (optim.crit == "WCC" && trwdth != trwdth.res)) {
      v <- switch(optim.crit,
                  RMS = RMS(a, ref, samp, B),
                  WCC = WCC(a, ref, samp, B, trwdth.res))
    }
  }

  ## calculate, or possibly re-calculate, quality of current solution
  if (try) {
    if (optim.crit == "WCC") {
      v <- WCC(a, ref, samp, B, trwdth.res)
    } else {
      if (smooth.param > 0) {
        samp <- t(apply(samp, 1, difsm, smooth.param))
        ref <- t(apply(ref, 1, difsm, smooth.param))
      }
      
      v <- RMS(a, ref, samp, B)
    }
  }

  ## back-transform coefficients
  w <- B %*% a
  a <- a/ncr^(0:(n-1))
  
  list(w = w, a = a, v = v)
} 
