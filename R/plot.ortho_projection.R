#' @title Plot method for an object of class \code{ortho_projection}
#' @description
#' Plots objects of class \code{ortho_projection}
#' @aliases plot.ortho_projection
#' @usage \method{plot}{ortho_projection}(x, col = "dodgerblue", ...)
#' @param x an object of class \code{ortho_projection} (as returned by \code{ortho_projection}).
#' @param col the color of the plots (default is "dodgerblue")
#' @param ... arguments to be passed to methods.
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @seealso \code{\link{ortho_projection}}
#' @importFrom graphics barplot
#' @export

plot.ortho_projection <- function(x, col = "dodgerblue", ...) {
  if (x$method == "pls") {
    x_variance <- x$variance$x_var
  } else {
    x_variance <- x$variance
  }

  if (x$pc_selection$method == "opc") {
    tpl <- x$opc_evaluation[, c(1, ncol(x$opc_evaluation))]
    if ("mean_standardized_rmsd_Yr" %in% colnames(tpl)) {
      ylab <- "mean of the standardized RMSD of all Y variables"
    }
    if (colnames(tpl)[2] %in% c("rmsd_Yr", "rmsd")) {
      ylab <- "RMSD of Yr"
    }
    if (colnames(tpl)[2] == "kappa") {
      ylab <- "kappa index"
    }
    plot(tpl,
      type = "p",
      ylab = ylab, pch = 1, col = col, ...
    )
    grid()
    segments(tpl[, 1], 0, tpl[, 1], tpl[, 2], col = col)
  } else {
    opar <- par("mfrow")
    on.exit(par(opar))
    
    o_mfrow <- par()$mfrow
    par(mfrow = c(1, 2))
    barplot(x_variance[grep("^explained_var", rownames(x_variance)), ],
      horiz = F,
      names.arg = colnames(x_variance), ylim = c(0, 1),
      ylab = "Explained variance", col = col, ...
    )
    barplot(x_variance[grep("cumulative", rownames(x_variance)), ],
      horiz = F,
      names.arg = colnames(x_variance), ylim = c(0, 1),
      ylab = "Explained variance (cumulative)", col = col, ...
    )
  }
}
