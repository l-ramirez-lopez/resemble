#' @title iterates over elements to print
#' @description internal. used for printing
#' @keywords internal
cat_iter <- function(x) {
  iter_x <- iter(x, by = "cell", recycle = TRUE)

  next_e <- function() {
    nextElem(iter_x)
  }
  obj <- next_e
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}

#' @title fits pls models at each r-local iteration
#' @description internal
#' @keywords internal
biter <- function(itersubs,
                  Xu,
                  Yu,
                  iter_sequence,
                  optimization,
                  ncomp,
                  tune,
                  p,
                  number,
                  scale,
                  max_iter,
                  tol,
                  allow_parallel) {
  "%mydo%" <- get("%do%")
  if (allow_parallel & getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }

  its <- NULL
  rdf <- foreach(
    idx = iter_sequence,
    its = itersubs
  ) %mydo% {
    # #
    # # Step 2 - sample a training data set of size k from K without replacement
    # #
    # selected.idx <- sample(k_idx,
    #                        size = k,
    #                        replace = FALSE)
    #
    # # retrieve the spectra and dependent variable from the SL for the training
    # set selection
    # X <- Xr[selected.idx,]
    # Y <- Yr[selected.idx]

    # dataset <- cbind(X, Y)

    #
    # Step 3 calibrate a PLS model using the selected SL samples
    #

    if (tune) {
      # perform a leave-group-out cross.validation of the selected SL (k)
      # observations to determine a suitable number of pls factors
      plsv <- pls_cv(
        x = its$x,
        y = its$y,
        ncomp = ncomp,
        method = "pls",
        center = TRUE,
        scale = scale,
        min_component = 1,
        p = p,
        number = number,
        group = its$group,
        retrieve = FALSE,
        tune = TRUE,
        max_iter = max_iter,
        tol = tol
      )

      selected_pls_c <- which.min(plsv$cv_results$rmse_cv)[1]
    } else {
      selected_pls_c <- ncomp
    }

    # build the model
    mod <- opls_get_basics(
      X = its$x,
      Y = its$y,
      ncomp = selected_pls_c,
      scale = scale,
      maxiter = max_iter,
      tol = tol
    )

    if (optimization == "response") {
      # rs <- Xu - (Xu %*% mod$projection_mat %*% mod$X_loadings)

      y_pred <- predict_opls(
        bo = mod$bo,
        b = mod$coefficients,
        ncomp = selected_pls_c,
        newdata = Xu,
        scale = scale,
        Xscale = mod$transf$Xscale
      )
      y_pred <- y_pred[, selected_pls_c]

      #
      # Step 4 - validate on the 'm' site specific observations
      #
      # y_pred <- as.numeric(predict(mod, Xu, ncomp=selected_pls_c))

      # It is possible we get NA predictions, catch this case and ignore this
      # observation test
      if (any(is.na(y_pred))) {

      } else {
        # rmse <- sqrt(sum((y_pred - Yu)^2)/(length(Yu) - 1))
        # rmse <- sqrt(sum((y_pred - Yu)^2)/(length(Yu)))
        rmse <- sqrt(mean((y_pred - Yu)^2))
      }
    }

    if (optimization == "reconstruction") {
      ## this must be implemented in Rcpp
      # xpred <- (Xu %*% mod$projection_mat) %*% mod$X_loadings
      # rmse <- sqrt(mean((Xu - xpred)^2))
      rmse <- reconstruction_error(Xu, mod$projection_mat, mod$X_loadings)[[1]]
    }

    # combine the validation results and the observations that were selected for this
    # model
    return(c(
      rmse_subset = rmse,
      sample = its$idx
    ))
  }

  # compile a data frame of all results of the sampling iterations
  rdf <- do.call("rbind", rdf)
  rmsesubset <- rdf[, "rmse_subset"]
  rdf <- rdf[, !colnames(rdf) %in% "rmse_subset"]
  return(list(
    rmsesubset = rmsesubset,
    sampleidx = rdf
  ))
}

#' @title iterator that subsets a matrix based on input indices for pls modeling
#' at each rs-local iteration
#' @description internal
#' @param x an input matrix of predictors
#' @param y a response variable
#' @param group a variable giving the groups/calsses to which the observations
#' belong to (used for avoiding pseudo-replication during validation)
#' @param indx the indices to retrieve
#' @return an object of \code{class} iterator giving a \code{list} with (a) a
#' subset of the input matrix with rows re-arranged according to `indx`
#' @details this is designed to be used within (parallel) loops to avoid sending
#' entire matrices to each core. It just sends the subset of that is required.
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
ithrssubsets <- function(x,
                         y,
                         group = NULL,
                         indx) {
  it_indx <- iter(indx, by = "column")

  if (is.null(group)) {
    nextEl <- function() {
      ss <- nextElem(it_indx)
      list(
        x = x[ss, , drop = FALSE],
        y = y[ss, , drop = FALSE],
        group = NULL,
        idx = ss
      )
    }
  } else {
    nextEl <- function() {
      ss <- nextElem(it_indx)
      list(
        x = x[ss, , drop = FALSE],
        y = y[ss, , drop = FALSE],
        group = factor(group[ss]),
        idx = ss
      )
    }
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}
