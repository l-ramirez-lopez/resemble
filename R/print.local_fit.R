#' @title Print method for an object of class \code{local_fit}
#' @description Prints the contents of an object of class \code{local_fit}
#' @usage \method{print}{local_fit}(x, ...)
#' @param x an object of class \code{local_fit}
#' @param ... not yet functional.
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#' @export

print.local_fit <- function(x, ...) {
  if (!x$method %in% c("pls", "wapls", "gpr")) {
    message("Method not recognized!")
  }

  if (x$method == "pls") {
    if (x$modified) {
      mss <- "Modified partial least squares (mpls)"
    } else {
      mss <- "Partial least squares (pls)"
    }
    cat(mss)
    cat("\nNumber of factors:", x$pls_c)
  }

  if (x$method == "wapls") {
    if (x$modified) {
      mss <- "Weighted average modified partial least squares (wampls)"
    } else {
      mss <- "Weighted average partial least squares (wapls)"
    }
    cat(mss)
    cat("\nMin. and max. number of factors: from", x$pls_c[["min_pls_c"]], "to", x$pls_c[["max_pls_c"]])
  }

  if (x$method == "gpr") {
    cat("Gaussian process with linear kernel/dot product (gpr)")
    cat("\nNoise:", x$noise_variance)
  }
}
