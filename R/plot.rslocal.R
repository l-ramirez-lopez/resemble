#' @title Plot method for an object of class \code{rslocal}
#' @description  Plots the content of an object of class \code{rslocal}
#' @aliases plot.rslocal
#' @usage \method{plot}{rslocal}(x, g = c("rmse", "n"), ...)
#' @param x an object of class \code{rslocal} (as returned by \code{rslocal}).
#' @param g a character vector indicating what results shall be plotted. 
#' Options are: \code{"rmse"} (for plotting the progress of the maximum RMSE found 
#' during the iterations) and/or \code{"n"} (for plotting the cumulative number 
#' of observations removed at each iteration).
#' @param ... some arguments to be passed to the plot methods.
#' @author Leonardo Ramirez-Lopez
#' @seealso \code{\link{rslocal}}
#' @export
###########################################################################

plot.rslocal <- function(x,
                         g = c("rmse", "n"), ...) {
 
  original_set <- par()$mfrow
  if (length(g) != 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
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
  } else{
    main <- "rs-local results"
  }

  grid_col <- rgb(0.3, 0.3, 0.3, 0.3)

  if ("rmse" %in% g) {
    do.call(plot, 
            c(list(x = x$iter_rmse$iteration, 
                   y = x$iter_rmse$max_rmse,
                   xlab = "# Iterations",
                   ylab = "Max. RMSE"),
              plot_dots))
    grid(col = grid_col, lty = 1)
  }

  if ("n" %in% g) {
    do.call(plot, 
            c(list(x = x$n_removed$iteration, 
                   y = x$n_removed$cummulative,
                   xlab = "# Iterations",
                   ylab = "Observations removed (cummulative)"),
              plot_dots))
    grid(col = grid_col, lty = 1)
   }
  mtext(main, outer = TRUE, cex = 2, line = -2)
  # par(ask = original_set)
  # op <- par(ask = original_set)
  on.exit(par(mfrow = original_set))
  # par(mfrow = pm)
  # title(main = "Memory-based learning results")
  # dev.flush()
}
