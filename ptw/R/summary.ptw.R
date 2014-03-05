summary.ptw <- function(object, ...)
{
  nsamp <- nrow(object$sample)
  nref <- nrow(object$reference)
  cat("PTW object:", object$warp.type,
      ifelse((object$warp.type == "individual" & nsamp > 1),
             "alignments of", "alignment of"), nsamp,
      ifelse(nsamp > 1, "samples on", "sample on"),
      nref, ifelse(nref > 1, "references.\n", "reference.\n"))
  cat("\nWarping coefficients:\n")
  print(object$warp.coef)
  cat("\nWarping criterion:", object$crit.type)
  cat(ifelse(object$warp.type == "individual" & nsamp > 1,
             "\nValues:", "\nValue:"), object$crit.value, "\n\n")
}

print.ptw <- function(x, ...)
{
  nsamp <- nrow(x$sample)
  nref <- nrow(x$reference)
  cat("PTW object:", x$warp.type,
      ifelse((x$warp.type == "individual" && nsamp > 1),
             "alignments of", "alignment of"), nsamp,
      ifelse(nsamp > 1, "samples on", "sample on"),
      nref, ifelse(nref > 1, "references.\n", "reference.\n"))
}
