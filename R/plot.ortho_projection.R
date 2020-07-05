#' @title Plot method for an object of class \code{ortho_projection}
#' @description
#'
#' \lifecycle{maturing}
#'
#' Plots the content pf an object of class \code{ortho_projection}
#' @aliases plot.ortho_projection
#' @usage \method{plot}{ortho_projection}(x, ...)
#' @param x an object of class \code{ortho_projection} (as returned by \code{ortho_projection}).
#' @param ... arguments to be passed to methods.
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @seealso \code{\link{ortho_projection}}
#' @importFrom graphics barplot
#' @export

plot.ortho_projection <- function(x, ...) {
  in.call <- match.call()$col
  if (is.null(in.call$col)) {
    col <- "dodgerblue"
  }

  if (x$method == "pls") {
    x_variance <- x$variance$x_var
  } else {
    x_variance <- x$variance
  }

  if (x$pc_selection$method == "opc") {
    tpl <- x$opc_evaluation[, c(1, ncol(x$opc_evaluation))]
    if ("mean_standardized_rmsd" %in% colnames(tpl)) {
      ylab <- "mean of the standardized RMSD of all Y variables"
    }
    if (colnames(tpl)[2] %in% c("rmsd_Yr", "rmsd")) {
      ylab <- "RMSD of Yr"
    }
    if (colnames(tpl)[2] == "kappa") {
      ylab <- "kappa index"
    }
    plot(tpl,
      type = "b",
      ylab = ylab, pch = 1, col = col, ...
    )
  }
  if (x$pc_selection$method == "cumvar") {
    barplot(x_variance[grep("cumulative", rownames(x_variance)), ],
      horiz = F,
      names.arg = colnames(x_variance), ylim = c(0, 1),
      ylab = "Explained variance (cummulative)", col = col, ...
    )
  }
  if (x$pc_selection$method %in% c("cumvar", "manual")) {
    x$variance
    barplot(x_variance[grep("^explained_var", rownames(x_variance)), ],
      horiz = F,
      names.arg = colnames(x_variance), ylim = c(0, 1),
      ylab = "Explained variance", col = col, ...
    )
  }
}
