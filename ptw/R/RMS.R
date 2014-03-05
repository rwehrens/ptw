RMS <- function(warp.coef, ref, samp, B)
{
  w <- B %*% warp.coef
  interp <- t(apply(samp, 1, function(x) interpol(w, x)))

  if (nrow(ref) == 1) ref <- c(ref)
  r <- interp - ref

  sqrt(mean(r^2, na.rm = TRUE))
}
