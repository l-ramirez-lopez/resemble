#' @title Plot method for an object of class \code{mbl}
#'
#' @description
#' Plots validation results and/or GH distance scores from an \code{mbl} object.
#'
#' @aliases plot.mbl
#'
#' @usage \method{plot}{mbl}(x, what = c("validation", "gh"),
#'     metric = "rmse", ncomp = c(1, 2), ...)
#'
#' @param x An object of class \code{mbl} (as returned by \code{\link{mbl}}).
#' @param what Character vector specifying what to plot. Options are
#'   \code{"validation"} (validation statistics) and/or \code{"gh"} (PLS scores
#'   used for GH distance computation). Default is both.
#' @param metric Character string specifying which validation statistic to plot.
#'   Options are \code{"rmse"}, \code{"st_rmse"}, or \code{"r2"}. Only used when
#'   \code{"validation"} is in \code{what}.
#' @param ncomp Integer vector of length 1 or 2 specifying which PLS components
#'   to plot. Default is \code{c(1, 2)}. Only used when \code{"gh"} is in
#'   \code{what}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @details
#' When plotting PLS scores (\code{what = "gh"}), the score matrix is
#' transformed from Euclidean to Mahalanobis space by multiplying by the
#' square root of the covariance matrix (computed via singular value
#' decomposition).
#'
#' @return Called for side effects (plotting). Returns \code{invisible(x)}.
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} and
#' Antoine Stevens
#'
#' @seealso \code{\link{mbl}}
#'
#' @examples
#' \dontrun{
#' # See mbl() examples for full workflow
#' 
#' library(prospectr)
#' data(NIRsoil)
#'
#' # Preprocess: detrend + first derivative with Savitzky-Golay
#' NIRsoil$spc_pr <- savitzkyGolay(
#'   (NIRsoil$spc),
#'   m = 1, p = 1, w = 7
#' ) |> standardNormalVariate()
#' NIRsoil$spc_pr <- sg_det
#' 
#' # Split data
#' test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso), ]
#' test_y <- NIRsoil$Ciso[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)]
#' train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso), ]
#' train_y <- NIRsoil$Ciso[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)]
#' 
#' mbl_result <- mbl(
#'   Xr = train_x,
#'   Yr = train_y,
#'   Xu = test_x,
#'   Yu = test_y,
#'   neighbors = neighbors_k(seq(40, 140, by = 20)),
#'   diss_method = diss_correlation(),
#'   fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15), 
#'   gh = TRUE
#' )
#' 
#' plot(mbl_result)
#' plot(mbl_result, what = "validation", metric = "r2")
#' plot(mbl_result, what = "gh", ncomp = c(2, 3))
#' }
#'
#' @export
plot.mbl <- function(x,
                     what = c("validation", "gh"),
                     metric = "rmse",
                     ncomp = c(1, 2),
                     ...) {
  # Validate arguments
  what <- match.arg(what, c("validation", "gh"), several.ok = TRUE)
  metric <- match.arg(metric, c("rmse", "st_rmse", "r2"))
  
  if (!inherits(x, "mbl")) {
    stop("'x' must be an object of class 'mbl'.", call. = FALSE)
  }
  
  # Save and restore graphics parameters
  
  opar <- par("mfrow", "mar")
  on.exit(par(opar))
  
  # Set up multi-panel if needed
  if (length(what) == 2 && !is.null(x$gh)) {
    par(mfrow = c(1, 2))
  }
  
  
  # Process plot arguments
  plot_args <- .process_plot_args(list(...))
  
  # Plot validation results
  if ("validation" %in% what) {
    .plot_validation(x, metric, plot_args)
  }
  
  # Plot GH scores
  if ("gh" %in% what) {
    if (is.null(x$gh)) {
      message("GH distance not available in this object.")
    } else {
      .plot_gh_scores(x, ncomp, plot_args)
    }
  }
  
  mtext(plot_args$main, outer = TRUE, cex = 2, line = -2)
  invisible(x)
}


# =============================================================================
# Internal plotting helpers
# =============================================================================

#' Process plot arguments with defaults
#' @keywords internal
.process_plot_args <- function(dots) {
  args <- dots
  
  # Extract or set main title
  if ("main" %in% names(args)) {
    main <- args$main
    args$main <- NULL
  } else {
    main <- "Memory-based learning results"
  }
  
  # Set defaults
  if (!"col.axis" %in% names(args)) {
    args$col.axis <- grey(0.3)
  }
  if (!"pch" %in% names(args)) {
    args$pch <- 16
  }
  
  # Remove arguments that are set internally
  internal_args <- c("col", "xlab", "ylab", "type", "ylim", "xlim")
  args <- args[!names(args) %in% internal_args]
  
  args$main <- main
  args
}


#' Plot validation results
#' @keywords internal
.plot_validation <- function(object, metric, plot_args) {
  # Collect validation results
  val_data <- list()
  colors <- character()
  
  if (!is.null(object$validation_results$nearest_neighbor_validation)) {
    val_data$NNv <- cbind(
      object$validation_results$nearest_neighbor_validation,
      val = "NNv"
    )
    colors <- c(colors, "dodgerblue")
  }
  
  if (!is.null(object$validation_results$local_cross_validation)) {
    val_data$local_cv <- cbind(
      object$validation_results$local_cross_validation,
      r2 = NA,
      val = "local_cv"
    )
    colors <- c(colors, "green4")
  }
  
  if (!is.null(object$validation_results$Yu_prediction_statistics)) {
    val_data$Yu <- cbind(
      object$validation_results$Yu_prediction_statistics,
      val = "Yu prediction"
    )
    colors <- c(colors, "red")
  }
  
  if (!is.null(object$validation_results$Yr_fitted_statistics)) {
    val_data$Yr <- cbind(
      object$validation_results$Yr_fitted_statistics,
      val = "Yr fitted"
    )
    colors <- c(colors, "orange")
  }
  
  if (length(val_data) == 0) {
    par(mfrow = c(1, 1))
    message("No validation results to plot.")
    return(invisible(NULL))
  }
  
  tpl <- do.call(rbind, val_data)
  col_names <- colnames(tpl)
  
  # Determine ID variable
  id_var <- if ("k" %in% col_names) "k" else "k_diss"
  
  # Select columns for metric
  keep_cols <- c(id_var, metric, "val")
  keep_cols <- keep_cols[keep_cols %in% col_names]
  tpl_subset <- data.frame(tpl)[, keep_cols, drop = FALSE]
  
  # Reshape to wide format
  to_plot <- reshape(
    tpl_subset,
    timevar = "val",
    idvar = id_var,
    direction = "wide"
  )
  colnames(to_plot) <- gsub("[.]", " ", colnames(to_plot))
  
  # Remove r2 for local_cv (not available)
  if (metric == "r2") {
    drop_col <- grepl("local_cv", colnames(to_plot))
    to_plot <- to_plot[, !drop_col, drop = FALSE]
    colors <- colors[colors != "green4"]
  }
  
  # Plot
  y_range <- range(to_plot[, -1], na.rm = TRUE)
  y_range[2] <- y_range[2] * 1.1
  
  do.call("matplot", c(
    list(
      x = to_plot[, 1],
      y = to_plot[, -1, drop = FALSE],
      type = "b",
      xlab = id_var,
      ylab = metric,
      ylim = y_range,
      col = colors
    ),
    plot_args[!names(plot_args) %in% "main"]
  ))
  
  grid(col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1, lwd = 1)
  mtext("Validation results", col = grey(0.3))
  
  legend(
    "topright",
    legend = colnames(to_plot)[-1],
    col = colors,
    pch = plot_args$pch,
    box.lty = 0,
    bg = NA
  )
  
  invisible(NULL)
}


#' Plot GH distance scores
#' @keywords internal
.plot_gh_scores <- function(object, ncomp, plot_args) {
  xr_scores <- object$gh$projection$scores
  n_components <- object$gh$projection$ncomp
  
  # Transform to Mahalanobis space if multicomponent
  if (n_components > 1) {
    xr_scores <- euclid_to_mahal(xr_scores)
    xr_scores <- sweep(xr_scores, 2, colMeans(xr_scores), "-")
  }
  
  # Split Xr and Xu scores
  xu_idx <- grep("Xu_", rownames(xr_scores))
  xr_idx <- grep("Xr_", rownames(xr_scores))
  
  xu_scores <- xr_scores[xu_idx, , drop = FALSE]
  xr_scores <- xr_scores[xr_idx, , drop = FALSE]
  
  # Colors
  xr_col <- rgb(0, 0, 0.4, 0.5)
  xu_col <- rgb(1, 0, 0, 0.5)
  
  if (n_components == 1) {
    .plot_gh_1d(xr_scores, xu_scores, xr_col, xu_col, plot_args)
  } else {
    .plot_gh_2d(xr_scores, xu_scores, ncomp, xr_col, xu_col, plot_args)
  }
  
  invisible(NULL)
}


#' Plot 1D GH scores
#' @keywords internal
.plot_gh_1d <- function(xr_scores, xu_scores, xr_col, xu_col, plot_args) {
  scores_combined <- c(xr_scores[, 1], xu_scores[, 1])
  set_labels <- c(rep("Xr", nrow(xr_scores)), rep("Xu", nrow(xu_scores)))
  
  df <- data.frame(
    index = seq_along(scores_combined),
    score = scores_combined,
    set = set_labels
  )
  df <- df[order(df$score), ]
  df$index <- seq_len(nrow(df))
  
  rng <- range(scores_combined)
  rng <- 1.2 * c(-max(abs(rng)), max(abs(rng)))
  
  do.call("plot", c(
    list(
      x = df$index[df$set == "Xr"],
      y = df$score[df$set == "Xr"],
      ylim = rng,
      col = xr_col,
      xlab = "Ordered PLS values",
      ylab = "PLS 1"
    ),
    plot_args[!names(plot_args) %in% "main"]
  ))
  
  points(
    df$index[df$set == "Xu"],
    df$score[df$set == "Xu"],
    col = xu_col,
    pch = plot_args$pch
  )
  
  mtext("Partial least squares scores", col = grey(0.3))
  grid(col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1, lwd = 1)
  legend("topleft", legend = c("Xr", "Xu"), col = c(xr_col, xu_col),
         pch = plot_args$pch, cex = 0.8, box.lty = 0, bg = NA)
}


#' Plot 2D GH scores
#' @keywords internal
.plot_gh_2d <- function(xr_scores, xu_scores, ncomp, xr_col, xu_col, plot_args) {
  rng <- 1.2 * range(xr_scores[, ncomp], xu_scores[, ncomp])
  
  xlab <- paste0("PLS ", ncomp[1], " (standardized)")
  ylab <- paste0("PLS ", ncomp[2], " (standardized)")
  
  do.call("plot", c(
    list(
      x = xr_scores[, ncomp[1]],
      y = xr_scores[, ncomp[2]],
      xlab = xlab,
      ylab = ylab,
      xlim = rng,
      ylim = rng,
      col = xr_col
    ),
    plot_args[!names(plot_args) %in% "main"]
  ))
  
  points(xu_scores[, ncomp, drop = FALSE], col = xu_col, pch = plot_args$pch)
  
  mtext("Partial least squares (Mahalanobis space)", col = grey(0.3))
  grid(col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1, lwd = 1)
  
  legend("topright", legend = c("Xr", "Xu"), col = c(xr_col, xu_col),
         pch = plot_args$pch, cex = 0.8, box.lty = 0, bg = NA)
  
  # Draw reference circles
  max_score <- ceiling(max(abs(xr_scores[, ncomp])))
  for (r in seq_len(max_score)) {
    theta <- seq(0, 2 * pi, length.out = 101)
    lines(r * cos(theta), r * sin(theta),
          col = rgb(0.3, 0.3, 0.3, 0.3), lty = 1, lwd = 0.5)
  }
}
