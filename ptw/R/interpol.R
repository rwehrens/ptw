interpol <- function (xout, y)
{                                                                             
  .C("R_approx", as.double(1:length(y)), as.double(y), as.integer(length(y)),
     xout = as.double(xout), as.integer(length(xout)), as.integer(1),
     as.double(NA), as.double(NA), as.double(0), NAOK = TRUE,
     PACKAGE = "ptw")$xout
}
