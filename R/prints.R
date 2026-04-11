# helpers
.divider <- function(width = 55L) {
  sys_width <- getOption("width")
  if (width > sys_width) width <- sys_width
  paste(rep("_", width), collapse = "")
}

.use_color <- function() {
  getOption("cli.num_colors", if (interactive()) 8L else 1L) > 1L
}

.col_blue <- function(x) {
  if (.use_color()) paste0("\033[34m", x, "\033[39m") else x
}

.col_bold_red <- function(x) {
  if (.use_color()) paste0("\033[1;31m", x, "\033[0m") else x
}

.truncate_call <- function(call_obj, width = 80L) {
  call_str <- deparse(call_obj, width.cutoff = width)
  if (length(call_str) > 3L) {
    call_str <- c(call_str[1:3], "...")
  }
  call_str
}

# -----------------------------------------------------------------------------
# print.liblex
# -----------------------------------------------------------------------------
#' @noRd
#' @export
print.liblex <- function(x, ...) {
  
  div <- .divider()
  
  has_coefficients <- !is.null(x$coefficients)
  
  # Header
  if (has_coefficients) {
    cat(.col_bold_red("--- liblex model library ---"), "\n")
    n_models <- length(x$coefficients$B0)
    cat(.col_blue("Models:    "), n_models, "\n")
    cat(.col_blue("Predictors:"), ncol(x$coefficients$B), "\n")
  } else {
    cat(.col_bold_red("--- liblex validation ---"), "\n")
    n_obs <- nrow(x$dissimilarity$dissimilarity)
    cat(.col_blue("Observations:"), n_obs, "\n")
  }
  cat(div, "\n")
  
  # Dissimilarity
  if (inherits(x$dissimilarity$diss_method, "diss_method")) {
    cat(.col_blue("Dissimilarity"), "\n")
    print(x$dissimilarity$diss_method)
    cat(div, "\n")
  }
  
  cat(.col_blue("Local fit method"), "\n")
  print(x$fit_method)
  cat(div, "\n")
  
  # Optimal parameters
  cat(.col_blue("Optimal parameters"), "\n")
  if (!is.null(x$optimal_params)) {
    if (inherits(x$neighbors, "neighbors_k"))
      cat("  k:    ", x$optimal_params$k, "\n")
    else
      cat("  diss threshold:    ", x$optimal_params$diss_threshold, "\n")
    cat("  ncomp:", x$optimal_params$ncomp["min"], "-", x$optimal_params$ncomp["max"], "\n")
  } else if (!is.null(x$best)) {
    cat("  k:    ", x$best$k, "\n")
    cat("  ncomp:", x$best$min_ncomp, "-", x$best$max_ncomp, "\n")
  }
  cat(div, "\n")
  
  # Validation statistics
  if (!is.null(x$best)) {
    cat(.col_blue("Nearest-neighbor validation"), "\n\n")
    val_df <- data.frame(
      rmse = round(x$best$rmse, 3),
      st_rmse = round(x$best$st_rmse, 3),
      me = round(x$best$me, 3),
      r2 = round(x$best$r2, 3)
    )
    print(val_df, row.names = FALSE)
    cat(div, "\n")
  }
  
  if (!has_coefficients) {
    cat("\nNote: Rebuild with mode = 'build' to enable predictions.\n")
  }
  
  invisible(x)
}

#' @noRd
#' @export
print.gesearch <- function(x, ...) {
  div <- .divider()
  val <- lapply(x$validation_results, function(xx) xx$results)
  
  # Header
  cat(.col_bold_red("--- gesearch results ---"), "\n")
  cat(.col_blue("Iterations:"), x$complete_iter, "\n")
  cat(.col_blue("Selected:  "), length(x$indices), "\n")
  cat(.col_blue("Removed:   "), max(x$n_removed$cumulative), "\n")
  cat(div, "\n")
  
  # Method info
  cat(.col_blue("Fit method"), "\n")
  print(x$fit_method)
  cat(div, "\n")
  
  # Training validation stats
  if (!is.null(val[[1L]]$train)) {
    cat(.col_blue("Validation (training set)"), "\n")
    for (i in seq_along(val)) {
      if (!is.null(names(val)) && names(val)[i] != "") {
        cat("Response:", names(val)[i], "\n")
      }
      print(as.data.frame(val[[i]]$train), digits = 3, row.names = FALSE)
      if (i < length(val)) cat("\n")
    }
    cat(div, "\n")
  }
  
  # Test validation stats
  if (!is.null(val[[1L]]$test)) {
    cat(.col_blue("Validation (test set)"), "\n")
    for (i in seq_along(val)) {
      if (!is.null(names(val)) && names(val)[i] != "") {
        cat("Response:", names(val)[i], "\n")
      }
      print(as.data.frame(val[[i]]$test), digits = 3, row.names = FALSE)
      if (i < length(val)) cat("\n")
    }
    cat(div, "\n")
  }
  
  invisible(x)
}

# -----------------------------------------------------------------------------
# print.mbl
# -----------------------------------------------------------------------------
#' @noRd
#' @export
print.mbl <- function(x, ...) {
  div <- .divider()
  val <- x$validation_results
  
  # Summary
  cat(.col_bold_red("--- mbl predictions ---"), "\n")
  cat(.col_blue("Predictions:"), x$n_predictions, "\n")
  cat(div, "\n")
  
  # Dissimilarity
  if (inherits(x$dissimilarities$diss_method, "diss_method")) {
    cat(.col_blue("Dissimilarity"), "\n")
    print(x$dissimilarities$diss_method)
    cat(div, "\n")
  }
  
  # Method info
  cat(.col_blue("Local fit method"), "\n")
  print(x$fit_method)
  cat(div, "\n")
  
  # Validation results
  if (!is.null(val$nearest_neighbor_validation)) {
    cat("\n")
    cat(.col_blue("Nearest-neighbor validation"), "\n\n")
    print(as.data.frame(val$nearest_neighbor_validation), digits = 3, row.names = FALSE)
    cat(div, "\n")
  }
  
  if (!is.null(val$local_cross_validation)) {
    cat(.col_blue("Local leave-group-out cross-validation"), "\n\n")
    print(as.data.frame(val$local_cross_validation), digits = 3, row.names = FALSE)
    cat(div, "\n")
  }
  
  if (!is.null(val$Yu_prediction_statistics)) {
    cat(.col_blue("Yu prediction statistics"), "\n\n")
    print(as.data.frame(val$Yu_prediction_statistics), digits = 3, row.names = FALSE)
    cat(div, "\n")
  }
  
  if (!is.null(val$Yr_fitted_statistics)) {
    cat(.col_blue("Yr fitted statistics"), "\n\n")
    print(as.data.frame(val$Yr_fitted_statistics), digits = 3, row.names = FALSE)
    cat(div, "\n")
  }
  
  invisible(x)
}

# =============================================================================
# Print method
# =============================================================================

#' @noRd
#' @export
#' @noRd
#' @export
print.resemble_model <- function(x, ...) {
  div <- .divider()
  
  # Header
  cat(.col_bold_red("--- Global resemble model ---"), "\n")
  cat(.col_blue("Method:      "), x$fit_method$fit_method, "\n")
  cat(.col_blue("Observations:"), x$n_obs, "\n")
  cat(.col_blue("Variables:   "), x$n_vars, "\n")
  cat(div, "\n")
  
  # Fit method details
  cat(.col_blue("Fit method"), "\n")
  print(x$fit_method)
  cat(div, "\n")
  
  # Cross-validation results
  if (!is.null(x$cv_results) && nrow(x$cv_results) > 0) {
    cat(.col_blue("Cross-validation"), "\n")
    if (inherits(x$fit_method, "fit_pls")) {
      best_idx <- which(x$cv_results$optimal)
      cat("  Best ncomp:", best_idx, "\n\n")
    }
    print(
      as.data.frame(x$cv_results[, names(x$cv_results) != "optimal"]),
      digits = 3,
      row.names = FALSE
    )
    cat(div, "\n")
  } else {
    cat(.col_blue("Validation:"), x$control$validation_type, "\n")
    cat(div, "\n")
  }
  
  invisible(x)
}
