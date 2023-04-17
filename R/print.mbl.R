#' @title Print method for an object of class \code{mbl}
#' @description Prints the content of an object of class \code{mbl}
#' @aliases print.mbl
#' @usage \method{print}{mbl}(x, ...)
#' @param x an object of class \code{mbl} (as returned by the \code{mbl} function).
#' @param ... arguments to be passed to methods (not functional).
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @keywords internal
#' @export

print.mbl <- function(x, ...) {
  val <- x$validation_results

  if (!is.null(val$nearest_neighbor_validation)) {
    nn_val_stats <- val$nearest_neighbor_validation
  } else {
    nn_val_stats <- NULL
  }
  if (!is.null(val$local_cross_validation)) {
    local_cv_stats <- val$local_cross_validation
  } else {
    local_cv_stats <- NULL
  }
  if (!is.null(val$Yu_prediction_statistics)) {
    yu_prediction_stats <- val$Yu_prediction_statistics
  } else {
    yu_prediction_stats <- NULL
  }
  if (!is.null(val$Yr_fitted_statistics)) {
    yr_fitted_statistics <- val$Yr_fitted_statistics
  } else {
    yr_fitted_statistics <- NULL
  }
  
  sys_width <- getOption("width")
  bar_width <- 55

  if (bar_width > sys_width) {
    bar_width <- sys_width
  }

  div <- paste(rep("_", bar_width), collapse = "")

  cat("\n")
  cat("Call:", "\n\n")
  print(x$call)
  cat("\n")
  cat(div, "\n")
  cat("\n", "Total number of observations predicted:", x$n_predictions, "\n")
  cat(div, "\n")

  if (!is.null(yr_fitted_statistics)) {
    cat("\n", "Statistics of the fitted Yr", "\n\n")
    print(yr_fitted_statistics, digits = 3)
    cat(div, "\n")
  }
  
  if (!is.null(nn_val_stats)) {
    cat("\n", "Nearest neighbor validation statistics", "\n\n")
    print(nn_val_stats, digits = 3)
    cat(div, "\n")
  }

  if (!is.null(local_cv_stats)) {
    cat("\n", "Average statistics of the local leave-group-out", "\n", "cross-validation", "\n\n")
    print(local_cv_stats, digits = 3)
    cat(div, "\n")
  }

  if (!is.null(yu_prediction_stats)) {
    cat("\n", "Statistics of the prediction of Yu", "\n\n")
    print(yu_prediction_stats, digits = 3)
  }
}
