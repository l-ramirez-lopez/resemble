#' @title Print method for an object of class \code{rslocal}
#' @description Prints the content of an object of class \code{rslocal}
#' @aliases print.rslocal
#' @usage \method{print}{rslocal}(x, ...)
#' @param x an object of class \code{rslocal} (as returned by the \code{rslocal} function).
#' @param ... arguments to be passed to methods (not functional).
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#' @export

print.rslocal <- function(x, ...) {
  val <- x$validation_results$results
  
  if (!is.null(val$train)) {
    train_stats <- val$train
  } else {
    train_stats <- NULL
  }
  if (!is.null(val$test)) {
    test_stats <- val$test
  } else {
    test_stats <- NULL
  }

  sys_width <- getOption("width")
  bar_width <- 55
  
  if (bar_width > sys_width) {
    bar_width <- sys_width
  }
  
  div <- paste(rep("_", bar_width), collapse = "")
    
  cat("\n")
  cat("Call:", "\n\n")
  print(attr(x, "call"))
  cat("\n")
  cat(div, "\n\n")
  cat("Iterations:", max(x$iter_rmse$iteration), "\n")
  cat("Total number of observations selected:", length(x$indices), "\n")
  cat("Total number of observations discarded:", max(x$n_removed$cummulative), "\n")
  cat("Final number of pls factors:", x$final_model$npls, "\n")
  
  cat(div, "\n")

  if (!is.null(train_stats)) {
    cat("\n", "Validation statistics for final training set", "\n\n")
    print(train_stats, digits = 3)
    cat(div, "\n")
  }

  if (!is.null(test_stats)) {
    cat("\n", "Validation statistics for test set", "\n\n")
    print(test_stats, digits = 3)
    cat(div, "\n")
  }
}
