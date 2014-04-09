RMS <- function(warp.coef, ref, samp, B, mode)
{
  w <- B %*% warp.coef

  if (mode == "backward") {
    interp <- t(apply(samp, 1, function(x) interpol(w, x)))
  } else {
    interp <- apply(samp,
                    1,
                    function(x) {
                      approx(w, x, xout = 1:length(x))$y
                    })
  }
  ## If samp is a one-row matrix interp will be a vector. Is this bad?

  if (nrow(ref) == 1) ref <- c(ref)
  r <- interp - ref

  sqrt(mean(r^2, na.rm = TRUE))
}
