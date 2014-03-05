ptw <- function (ref, samp, selected.traces,
                 init.coef = c(0, 1, 0), try = FALSE,
		 warp.type = c("individual", "global"),
		 optim.crit = c("WCC", "RMS"),
		 smooth.param = ifelse(try, 0, 1e05),
		 trwdth = 20, trwdth.res = trwdth,
                 verbose = FALSE,
                 ... )
{
  optim.crit <- match.arg(optim.crit)
  warp.type <- match.arg(warp.type)

  if (is.vector(ref)) ref <- matrix(ref, nrow = 1)
  if (is.vector(samp)) samp <- matrix(samp, nrow = 1)
  if (nrow(ref) > 1 && nrow(ref) != nrow(samp))
    stop("The number of references does not equal the number of samples")
  if (nrow(samp) == 1) warp.type <- "individual"
  
  r <- nrow(samp)
  
  if (!missing(selected.traces)) {
    samp <- samp[selected.traces,, drop=FALSE]  

    if (nrow(ref) > 1)
      ref <- ref[selected.traces,, drop=FALSE]
          
  }

  if (is.vector(init.coef)) init.coef <- matrix(init.coef, nrow = 1)
  if (warp.type == "global") {
    if (nrow(init.coef) != 1)
      stop("Only one warping function is allowed with global alignment.")
  } else {
    if (nrow(init.coef) != nrow(samp))
      if (nrow(init.coef) == 1) {
        init.coef <- matrix(init.coef, byrow = TRUE,
                         nrow = nrow(samp), ncol = length(init.coef))
      } else {
        stop("The number of warping functions does not match the number of samples")
      }
  }
  
  if (warp.type == "individual") {
    w <- matrix(0, nrow(samp), ncol(ref))   
    a <- matrix(0, nrow(samp), ncol(init.coef))
    v <- rep(0, nrow(samp))
    warped.sample <- matrix(NA, nrow=nrow(samp), ncol=ncol(samp))
    
    for (i in 1:nrow(samp)) {
      if (verbose & nrow(samp) > 1)
        cat(ifelse(nrow(ref) == 1,
                   paste("Warping sample", i, "with the reference \n"),
                   paste("Warping sample", i, "with reference \n", i)))

      if (nrow(ref) == 1) {
        rfrnc <- ref
      } else {
        rfrnc <- ref[i, , drop = FALSE]
      }
      switch(optim.crit,
             RMS = {
               quad.res <- pmwarp(rfrnc,
                                  samp[i, , drop = FALSE],
                                  optim.crit, init.coef[i,], try = try,
                                  smooth.param = smooth.param, ...)
             },
             WCC = {
               quad.res <- pmwarp(rfrnc,
                                  samp[i, , drop = FALSE],
                                  optim.crit, init.coef[i,], try = try,
                                  trwdth = trwdth, trwdth.res = trwdth.res,
                                  ...)
             })
      
      w[i, ] <- quad.res$w
      a[i, ] <- quad.res$a
      v[i] <- quad.res$v
      warped.sample[i, ] <- interpol(w[i, ], samp[i, ])
    }
  } else {
    warped.sample <- matrix(NA, nrow=nrow(samp), ncol=ncol(samp))
    
    if (nrow(ref)==1) 
      ref <- matrix(ref, nrow = nrow(samp), 
      			ncol = ncol(ref), byrow = TRUE)
    
    if (verbose) {
      if (nrow(ref) == 1) {
        cat("Simultaneous warping of samples with reference... \n")
      } else {
        cat("Simultaneous warping of samples with references... \n")
      }
    }

    switch(optim.crit,
           RMS = {
             quad.res <- pmwarp(ref, samp, optim.crit, c(init.coef), try = try,
                                smooth.param = smooth.param,
                                ...)
           },
           WCC = {
             quad.res <- pmwarp(ref, samp, optim.crit, c(init.coef), try = try,
                                trwdth = trwdth, trwdth.res = trwdth.res,
                                ...)
           })
    
    w <- t(as.matrix(quad.res$w))
    a <- t(as.matrix(quad.res$a))
    v <- quad.res$v
    warped.sample <- t(sapply(1:nrow(samp),
                         function(i) {
                           interpol(w, samp[i,])
                         }))
  }
    
  if (verbose) cat("\nFinished.\n")  

  result <-list(reference = ref, sample = samp,
                warped.sample = warped.sample,
                warp.coef = a, warp.fun = w,
                crit.value = v, optim.crit = optim.crit,
                warp.type = warp.type)
  class(result) <- "ptw"
  result
}
