#' @title Plot method for ortho_projection objects
#' 
#' @description
#' Plots variance explained or OPC evaluation results for objects of class
#' \code{ortho_projection}.
#'
#' @param x An object of class \code{ortho_projection} (as returned by
#'   \code{\link{ortho_projection}}).
#' @param col Color for the plot elements. Default is \code{"#3B82F6"}.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @details
#' The plot type depends on the component selection method used:
#' \itemize{
#'   \item For \code{\link{ncomp_by_opc}}: displays RMSD (or kappa) as a 
#'     function of the number of components.
#'   \item For other methods: displays individual and cumulative explained 
#'     variance in a two-panel layout.
#' }
#'
#' @return Invisible \code{NULL}. Called for its side effect of producing a plot.
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} and
#' Antoine Stevens
#'
#' @seealso \code{\link{ortho_projection}}, \code{\link{ncomp_by_opc}}
#'
#' @importFrom graphics barplot grid par plot segments
#' @export
plot.ortho_projection <- function(x, col = "#3B82F6", ...) {
  
  if (inherits(x$ncomp_method, "ncomp_by_opc")) {
    .plot_opc_evaluation(x, col = col, ...)
  } else {
    .plot_variance(x, col = col, ...)
  }
  
  invisible(NULL)
}

# --- internal: OPC evaluation plot --------------------------------------------
.plot_opc_evaluation <- function(x, col, ...) {
  opc_eval <- x$opc_evaluation
  tpl <- opc_eval[, c(1, ncol(opc_eval))]
  metric <- colnames(tpl)[2]
  
  ylab <- switch(
    metric,
    mean_standardized_rmsd_Yr = "Mean standardized RMSD",
    rmsd_Yr = "RMSD",
    rmsd = "RMSD",
    kappa = "Kappa",
    metric
  )
  
  plot(
    tpl,
    type = "p",
    pch = 16,
    col = col,
    ylab = ylab,
    xlab = "Number of components",
    las = 1,
    ...
  )
  grid(col = "#33415580", lty = 1)
  segments(tpl[, 1], 0, tpl[, 1], tpl[, 2], col = col)
  
  # Mark optimal
  opt_idx <- which.min(tpl[, 2])
  if (metric == "kappa") {
    opt_idx <- which.max(tpl[, 2])
  }
  points(tpl[opt_idx, 1], tpl[opt_idx, 2], pch = 16, col = "#F59E0B", cex = 1.5)
}

# --- internal: variance plot --------------------------------------------------
.plot_variance <- function(x, col, ...) {
  x_var <- x$variance$x_var
  
  opar <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit(par(opar))
  
  # Individual variance
  ind_var <- x_var[grep("^explained_var", rownames(x_var)), ]
  barplot(
    ind_var,
    horiz = FALSE,
    names.arg = colnames(x_var),
    ylim = c(0, 1),
    ylab = "Explained variance",
    xlab = "Component",
    col = col,
    border = NA,
    las = 1,
    ...
  )
  
  # Cumulative variance
  cum_var <- x_var[grep("cumulative", rownames(x_var)), ]
  barplot(
    cum_var,
    horiz = FALSE,
    names.arg = colnames(x_var),
    ylim = c(0, 1),
    ylab = "Cumulative explained variance",
    xlab = "Component",
    col = col,
    border = NA,
    las = 1,
    ...
  )
}
