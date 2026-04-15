#' @title Plot method for an object of class \code{gesearch}
#' @description Plots objects of class \code{gesearch}
#' @aliases plot.gesearch
#' @usage \method{plot}{gesearch}(x, g = c("rmse", "n"), ...)
#' @param x an object of class \code{gesearch} (as returned by \code{gesearch}).
#' @param g a character vector indicating what results shall be plotted.
#' Options are: \code{"rmse"} (for plotting the progress of the maximum RMSE 
#' found during the iterations) or \code{"n"} (for plotting the cumulative 
#' number of observations removed at each iteration).
#' @param ... some arguments to be passed to the plot methods.
#' @author Leonardo Ramirez-Lopez
#' @seealso \code{\link{gesearch}}
#' @export
###########################################################################

plot.gesearch <- function(x, g = c("rmse", "n"), ...) {
  
  if (length(g) > 1 | !any(g %in% c("rmse", "n"))) {
    stop("g can only be either 'rmse' or 'n'")
  }
  
  opar <- par("mfrow", "mar")
  on.exit(par(opar))
  
  if (g == "rmse") {
    nbox <- ceiling(sqrt(ncol(x$iter_weakness$maximum) - 1))
    op <- par(mfrow = c(nbox, nbox), mar = c(4, 5.5, 2.5, 2))
  }
  
  plot_dots <- list(...)
  
  if (!"col" %in% names(plot_dots)) {
    plot_dots$col <- "dodgerblue"
  }
  
  if (!"type" %in% names(plot_dots)) {
    plot_dots$type <- "b"
  }
  
  if ("xlab" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "xlab"]
  }
  
  if ("ylab" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "ylab"]
  }
  
  if ("main" %in% names(plot_dots)) {
    main <- plot_dots$main
    plot_dots <- plot_dots[!names(plot_dots) %in% "main"]
  } else {
    main <- "gesearch results"
  }
  
  grid_col <- rgb(0.3, 0.3, 0.3, 0.3)
  
  if ("rmse" %in% g) {
    
    pnames <- colnames(x$iter_weakness$maximum)[-1]
    for (i in 1:length(pnames)) {
      ith_ylab <- gsub("_", " ", gsub("rmse_", "", pnames[i]))
      ith_ylab <- paste0("Max. RMSE\n", ith_ylab)
      ith_col <- 1 + i
      do.call(
        plot,
        c(
          list(
            x = x$iter_weakness$maximum[[1]],
            y = x$iter_weakness$maximum[[ith_col]],
            xlab = "# Iterations",
            ylab = ith_ylab
          ),
          plot_dots
        )
      )
      grid(col = grid_col, lty = 1)
    }
  }
  
  if ("n" %in% g) {
    do.call(
      plot,
      c(
        list(
          x = x$n_removed$iteration,
          y = x$n_removed$cummulative,
          xlab = "# Iterations",
          ylab = "Observations removed (cummulative)"
        ),
        plot_dots
      )
    )
    grid(col = grid_col, lty = 1)
  }
  mtext(main, outer = TRUE, cex = 2, line = -2)
  # par(ask = original_set)
  # op <- par(ask = original_set)
  # on.exit(par(mfrow = original_set))
  # par(mfrow = pm)
  # title(main = "Memory-based learning results")
  # dev.flush()
}
