#' @title Plot method for an object of class \code{mbl}
#' @description
#' Plots the content of an object of class \code{mbl}
#' @aliases plot.mbl
#' @usage \method{plot}{mbl}(x, g = c("validation", "gh"), param = "rmse", pls_c = c(1,2), ...)
#' @param x an object of class \code{mbl} (as returned by \code{mbl}).
#' @param g a character vector indicating what results shall be plotted.
#' Options are: \code{"validation"} (for plotting the validation results) and/or
#' \code{"gh"} (for plotting the pls scores used to compute the GH distance.
#' See details).
#' @param param a character string indicating what validation statistics shall be
#' plotted. The following options are available: \code{"rmse"}, \code{"st_rmse"}
#' or \code{"r2"}. These options only available if the \code{mbl} object contains
#' validation results.
#' @param pls_c a numeric vector of length one or two indicating the pls factors to be
#' plotted. Default is \code{c(1, 2)}. It is only available if \code{"gh"} is
#' specified in the \code{g} argument.
#' @param ... some arguments to be passed to the plot methods.
#' @details
#' For plotting the pls scores from the pls score matrix (of more than one column),
#' this matrix is first transformed from the Euclidean space to the Mahalanobis
#' space. This is done by multiplying the score matrix by the root square of
#' its covariance matrix. The root square of this matrix is estimated using a
#' singular value decomposition.
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @seealso \code{\link{mbl}}
#' @examples
#' \donttest{
#' library(prospectr)
#'
#' data(NIRsoil)
#'
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' Xu <- Xu[!is.na(Yu), ]
#' Yu <- Yu[!is.na(Yu)]
#'
#' Xr <- Xr[!is.na(Yr), ]
#' Yr <- Yr[!is.na(Yr)]
#'
#' ctrl <- mbl_control(validation_type = "NNv")
#'
#' ex_1 <- mbl(
#'   Yr = Yr, Xr = Xr, Xu = Xu,
#'   diss_method = "cor",
#'   diss_usage = "none",
#'   gh = TRUE,
#'   mblCtrl = ctrl,
#'   k = seq(50, 250, 30)
#' )
#'
#' plot(ex_1)
#' plot(ex_1, g = "gh", pls_c = c(2, 3))
#' }
#' @export
###########################################################################
## History:
## 23.04.2014 Leo     Plot function when the data is not centred now
##                    draws the circles around the actual centre
## 12.04.2015         When the circle was plotted there was an small
##                    gap in it. This was fixed by modifiying the pntCirc
##                    function
## 19.08.2016         10.08.2016 A bug in plot.mbl was corrected.
##                    It was not possible to plot mbl results when the
##                    k.diss argument (threshold distances) was used in
##                    the mbl function.
## 22.06.2020         Modified to match the new otputs of mbl
##                    option "pca" was replaced by option "gh" which plots the
##                    pls projection used for computing the gh distance in mbl
plot.mbl <- function(x,
                     g = c("validation", "gh"),
                     param = "rmse",
                     pls_c = c(1, 2), ...) {
  opar <- par("mfrow", "mar")
  on.exit(par(opar))

  if (length(g) != 1 & !is.null(x$gh)) {
    op <- par(mfrow = c(1, 2))
  }

  plot_dots <- list(...)
  if ("main" %in% names(plot_dots)) {
    main <- plot_dots$main
    plot_dots <- plot_dots[!names(plot_dots) %in% "main"]
  } else {
    main <- "Memory-based learning results"
  }

  if (!"col.axis" %in% names(plot_dots)) {
    plot_dots$col.axis <- grey(0.3)
  }

  if (!"pch" %in% names(plot_dots)) {
    plot_dots$pch <- 16
  }

  if ("col" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "col"]
  }

  if ("xlab" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "xlab"]
  }

  if ("ylab" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "ylab"]
  }

  if ("type" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "type"]
  }

  if ("ylim" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "ylim"]
  }

  if ("xlim" %in% names(plot_dots)) {
    plot_dots <- plot_dots[!names(plot_dots) %in% "xlim"]
  }

  object <- x
  # pm <- par()$mfrow

  if ("validation" %in% g) {
    col <- NULL
    if (!is.null(object$validation_results$nearest_neighbor_validation)) {
      nn_val_stats <- cbind(object$validation_results$nearest_neighbor_validation, val = "NNv")
      col <- c(col, "dodgerblue")
    } else {
      nn_val_stats <- NULL
    }
    if (!is.null(object$validation_results$local_cross_validation)) {
      local_cv_stats <- cbind(object$validation_results$local_cross_validation, r2 = NA, val = "local_cv")
      col <- c(col, "green4")
    } else {
      local_cv_stats <- NULL
    }
    if (!is.null(object$validation_results$Yu_prediction_statistics)) {
      yu_prediction_stats <- cbind(object$validation_results$Yu_prediction_statistics, val = "Yu prediction")
      col <- c(col, "red")
    } else {
      yu_prediction_stats <- NULL
    }

    tpl <- rbind(nn_val_stats, local_cv_stats, yu_prediction_stats)

    if (is.null(tpl)) {
      par(mfrow = c(1, 1))
      message("No validation results to plot")
    } else {
      # par(mfrow = c(1, length(g)))
      dtn <- colnames(tpl)
      opt <- c("rmse", "st_rmse", "r2")
      dt <- !is.element(dtn, opt[!is.element(opt, param)])
      idv <- ifelse("k" %in% colnames(tpl), "k", "k_diss")
      dt <- as.logical(dt * (!dtn %in% "p_bounded"))

      to_plot <- data.frame(tpl)[, dt] %>%
        stats::reshape(
          timevar = "val",
          idvar = idv,
          direction = "wide"
        )
      colnames(to_plot) <- gsub("[.]", " ", colnames(to_plot))

      if (param == "r2") {
        to_plot <- to_plot[, !colnames(to_plot) == "r2_local_cv"]
        col <- col[!col == "green4"]
      }
      do.call("matplot", c(
        list(
          x = to_plot[, 1],
          y = to_plot[, -1],
          type = "b",
          xlab = dtn[1],
          ylab = param,
          ylim = c(min(to_plot[, -1]), 1.1 * max(to_plot[, -1])),
          col = col
        ),
        plot_dots
      ))
      grid(
        nx = NULL, ny = NULL, col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1,
        lwd = 1, equilogs = TRUE
      )
      mtext("Validation results", col = grey(0.3))
      # Adding a legend
      legend("topright",
        legend = colnames(to_plot[, -1, drop = FALSE]), bg = NA,
        pch = plot_dots$pch,
        col = col, box.lty = 0
      )
      # is.element("ggplot2", installed.packages()[,"Package"])
    }
  }

  if ("gh" %in% g & !is.null(x$gh)) {
    xr_scores <- object$gh$projection$scores
    if (object$gh$projection$n_components > 1) {
      xr_scores <- euclid_to_mahal(xr_scores)
      mean_xr_scores <- colMeans(xr_scores)
      xr_scores <- sweep(xr_scores, MARGIN = 2, FUN = "-", STATS = mean_xr_scores)
    }

    xu_scores <- xr_scores[grep("Xu_", rownames(xr_scores)), , drop = FALSE]
    xr_scores <- xr_scores[grep("Xr_", rownames(xr_scores)), , drop = FALSE]
    xr_col <- rgb(0, 0, 0.4, 0.5)
    xu_col <- rgb(1, 0, 0, 0.5)
    circle_col <- rgb(0, 0.6, 0.6, 0.6)

    if (object$gh$projection$n_components == 1) {
      rng <- range(xr_scores, xu_scores)
      rng <- 1.2 * c(-max(abs(rng)), max(abs(rng)))
      tp <- c(xr_scores[, 1], xu_scores[, 1])
      tp <- data.frame(
        index = 1:length(tp), tp = tp,
        set = c(
          rep("Xr", nrow(xr_scores)),
          rep("Xu", nrow(xu_scores))
        )
      )
      tp <- tp[order(tp$tp), ]
      tp$index <- 1:length(tp$index)

      do.call("plot", c(
        list(
          x = tp[tp$set == "Xr", 1],
          y = tp[tp$set == "Xr", 2],
          ylim = rng,
          col = xr_col,
          ylab = "pls 1", xlab = "Ordered pls values"
        ),
        plot_dots
      ))

      mtext("Partial least squares scores", col = grey(0.3))
      points(tp[tp$set == "Xu", 1:2], xlim = rng, ylim = rng, col = xu_col, pch = plot_dots$pch)
      grid(
        nx = NULL, ny = NULL, col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1,
        lwd = 1, equilogs = TRUE
      )
      legend("topleft",
        legend = c("Xr", "Xu"),
        col = c(xr_col, xu_col), pch = 16, cex = 0.8, box.lty = 0, bg = NA
      )
    } else {
      rng <- 1.2 * range(xr_scores[, pls_c], xu_scores[, pls_c])

      xl <- paste0("pls ", pls_c[1], " (standardized)")
      yl <- paste0("pls ", pls_c[2], " (standardized)")

      do.call("plot", c(
        list(
          x = xr_scores[, pls_c[1]],
          y = xr_scores[, pls_c[2]],
          xlab = xl,
          ylab = yl,
          xlim = rng,
          ylim = rng,
          col = xr_col
        ),
        plot_dots
      ))

      mtext("Partial least squares (Mahalanobis space)", col = grey(0.3))


      points(xu_scores[, pls_c, drop = FALSE], col = xu_col, pch = plot_dots$pch)
      grid(
        nx = NULL, ny = NULL, col = rgb(0.3, 0.3, 0.3, 0.1), lty = 1,
        lwd = 1, equilogs = TRUE
      )
      legend("topright",
        legend = c("Xr", "Xu"),
        col = c(xr_col, xu_col), pch = plot_dots$pch, cex = 0.8, box.lty = 0,
        bg = NA
      )

      pntCirc <- function(r) {
        n <- 100
        a <- matrix(0, n + 1, 2)
        for (i in 1:n) {
          pnts <- (c(cos(2 * pi / n * i) * r, sin(2 * pi / n * i) * r))
          a[i, ] <- pnts
        }
        a[i + 1, ] <- a[1, ]
        return(a)
      }

      for (i in 1:ceiling(max(xr_scores[, pls_c]))) {
        crc <- pntCirc(i)
        crc <- rbind(crc, crc[1, ])
        lines(crc, col = rgb(0.3, 0.3, 0.3, 0.3), lty = 1, lwd = 0.5)
      }
    }
  } else {
    message("gh not available in this object")
  }
  mtext(main, outer = TRUE, cex = 2, line = -2)
  # par(ask = original_set)
  # op <- par(ask = original_set)
  # par(mfrow = pm)
  # title(main = "Memory-based learning results")
  # dev.flush()
}
