#' @title Print method for an object of class \code{gesearch}
#' @description Prints the content of an object of class \code{gesearch}
#' @aliases print.gesearch
#' @usage \method{print}{gesearch}(x, ...)
#' @param x an object of class \code{gesearch} (as returned by the \code{gesearch} function).
#' @param ... arguments to be passed to methods (not functional).
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#' @export

print.gesearch <- function(x, ...) {
  val <- lapply(x$validation_results, FUN = function(xx) xx$results)
  
  if (!is.null(val[[1]]$train)) {
    train_val <- lapply(val, FUN = function(xx) xx$train[, -1])
    # train_val <- do.call("cbind", train_val)
    # colnames(train_val) <- gsub("\\.", "_", colnames(train_val))
    # train_stats <- cbind(val[[1]]$train[,1], train_val)
  } else {
    train_val <- NULL
  }
  if (!is.null(val[[1]]$test)) {
    test_val <- lapply(val, FUN = function(xx) xx$test[, -1])
    # test_val <- do.call("cbind", test_val)
    # colnames(test_val) <- gsub("\\.", "_", colnames(test_val))
    # test_stats <- cbind(val[[1]]$test[,1], test_val)
  } else {
    test_val <- NULL
  }
  
  sys_width <- getOption("width")
  bar_width <- 55
  
  if (bar_width > sys_width) {
    bar_width <- sys_width
  }
  
  div <- paste(rep("_", bar_width), collapse = "")
  
  if (!is.null(names(x$final_models))) {
    response_names <- paste0("(", paste0(names(x$final_models), collapse = ", "), ")")
  } else {
    response_names <- NULL
  }
  
  cat("\n")
  cat("Call:", "\n\n")
  print(attr(x, "call"))
  cat("\n")
  cat(div, "\n\n")
  cat("Iterations:", max(x$complete_iter), "\n")
  cat("Number of response variables:", length(x$final_models), response_names, "\n")
  cat("Total number of observations selected:", length(x$indices), "\n")
  cat("Total number of observations discarded:", max(x$n_removed$cummulative), "\n")
  cat("Final number of pls factors:", x$final_model$npls, "\n")
  
  cat(div, "\n")
  
  if (!is.null(train_val)) {
    cat("\n", "--- Validation statistics for final training set ---", "\n\n")
    tr <- sapply(1:length(train_val), FUN = function(x, i) {
      if (!is.null(names(x))) {
        cat("Response:",  names(x)[i], "\n")
        print(x[[i]], digits = 3)
        cat("\n")
      } else {
        print(x[[i]], digits = 3)
      }
      
    },
    x = train_val)
    cat(div, "\n")
  }
  
  if (!is.null(test_val)) {
    cat("\n", "--- Validation statistics for test set ---", "\n\n")
    vl <- sapply(1:length(test_val), FUN = function(x, i) {
      if (!is.null(names(x))) {
        cat("Response:",  names(x)[i], "\n")
        print(x[[i]], digits = 3)
        cat("\n")
      } else {
        print(x[[i]], digits = 3)
      }
    },
    x = test_val)
    cat(div, "\n")
  }
}
