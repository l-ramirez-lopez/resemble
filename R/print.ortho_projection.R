#' @title Print method for an object of class \code{ortho_projection}
#' @description Prints the contents of an object of class \code{ortho_projection}
#' @aliases print.ortho_projection
#' @usage \method{print}{ortho_projection}(x, ...)
#' @param x an object of class \code{ortho_projection} (as returned by the
#' \code{ortho_projection} function).
#' @param ... arguments to be passed to methods (not yet functional).
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#' @export

print.ortho_projection <- function(x, ...) {
  cat("\n", "Method: ", x$method)
  cat("\n", "Number of components retained: ", x$n_components, "\n")
  cat(" Number of observations and number of original variables: ", c(nrow(x$scores), ncol(x$X_loadings)), "\n")
  cat(paste0(" Original variance in Xr: ",  round(x$variance$original_x_var, 3), "\n"))
  
  if (any(grepl("Xu", rownames(x$scores)))) {
    mss <- c("\n", "Explained variance in X {Xr; Xu}: \n")
  } else {
    mss <- c("\n", "Explained variance in Xr: \n")
  }
  
  cat("\n", "Explained variances, ratio of explained variance, cumulative explained variance:", "\n")
  cat(mss)
  print(x$variance$x_var, digits = 3)
  
  if (x$method %in% c("pls", "mpls")) {
    cat("\n", "Explained variance in Yr: \n")
    print(x$variance$y_var, digits = 3)
  }
}
