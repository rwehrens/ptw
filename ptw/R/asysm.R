asysm <- function (y, lambda = 1e+07, p = 0.001, eps = 1e-8) {
  ny <- length(y)
  w <- rep(1, ny)
  z <- d <- c <- e <- rep(0, length(y))
  R <- TRUE
  it <- 0
  while (R) {
    z <- .C("smooth2", as.double(w), as.double(y), as.double(z),
            as.double(lambda), as.integer(length(y)), as.double(d),
            as.double(c), as.double(e), PACKAGE = "ptw")[[3]]
    w0 <- w
    w <- p * (y - z > eps) + (1 - p) * (y - z <= eps)
    R <- sum(abs(w - w0)) > 0
  }
  z
} 