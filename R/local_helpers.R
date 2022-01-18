#' @title A function to get the neighbor information
#' @description This fucntion gathers information of all neighborhoods of the
#' \code{Xu} observations found in \code{Xr}. This information is equired during
#' local regressions.
#' @details
#' For local pca and pls distances, the local dissimilarity matrices are not
#' computed as it is cheaer to compute them during the local regressions.
#' Instead the global distances (required for later local dissimilarity matrix
#' computation are output)
#' @keywords internal
get_neighbor_info <- function(Xr, Xu, diss_method, Yr = NULL,
                              k = NULL, k_diss = NULL, k_range = NULL,
                              spike = NULL, pc_selection,
                              return_dissimilarity,
                              center, scale, gh,
                              diss_usage, allow_parallel = FALSE,
                              ...) {
  ortho_diss_methods <- c("pca", "pca.nipals", "pls")
  k_max <- NULL
  if (!is.null(k)) {
    k <- sort(k)
    k_max <- max(k)
  }

  if (!is.null(k_diss)) {
    k_diss <- sort(k_diss)
    k_diss_max <- max(k_diss)
    # k_min_range <- min(k_range)
    # k_max_range <- max(k_range)
  } else {
    k_diss_max <- k_range <- NULL
  }

  neighbor_dots <- list(...)
  # local dissimilarities option is always
  # deactivated, as this will be calculated inside the local loops for
  # modeling. This Saves memory
  is_local <- isTRUE(neighbor_dots$.local)
  if (is_local & diss_method %in% ortho_diss_methods) {
    neighbor_dots$.local <- FALSE
    if (is.null(neighbor_dots$pre_k)) {
      stop("'pre_k' argument is missing, it is required when .local = TRUE")
    }
    k_max <- neighbor_dots$pre_k
    k_diss_max <- k_range <- k_diss <- NULL
  }

  ifelse(is_local, return_diss <- FALSE, return_diss <- return_dissimilarity)
  # use do.call to be able to change local to FALSE in ... (if provided)
  info_neighbors <- do.call(
    search_neighbors,
    c(
      list(
        Xr = Xr, Xu = Xu, diss_method = diss_method,
        Yr = Yr,
        k = k_max,
        k_diss = k_diss_max, k_range = k_range,
        spike = spike,
        pc_selection = pc_selection,
        return_projection = TRUE,
        return_dissimilarity = return_diss,
        center = center, scale = scale,
        gh = gh
      ),
      neighbor_dots
    )
  )

  # if (!is.null(k_diss)) {
  #   neighbors_diss <- neighbors <- NULL
  #   neighbors[[length(k_diss)]] <- info_neighbors$neighbors
  #   neighbors_diss[[length(k_diss)]] <- info_neighbors$neighbors_diss
  #
  #   for (i in (length(k_diss) - 1):1) {
  #     neighbors[[i]] <- info_neighbors$neighbors
  #     neighbors_diss[[i]] <- info_neighbors$neighbors_diss
  #     i_is_neigh <- info_neighbors$neighbors_diss <= k_diss[i]
  #     i_n_neigh <- colSums(i_is_neigh, na.rm = TRUE)
  #     i_n_neigh[i_n_neigh < k_min_range] <- k_min_range
  #     i_n_neigh[i_n_neigh > k_max_range] <- k_max_range
  #
  #     list_neigh <- Map(':', rep(1, length(i_n_neigh)), i_n_neigh)
  #     sqnc <- seq(1, nrow(i_is_neigh) * ncol(i_is_neigh), by = nrow(i_is_neigh)) - 1
  #     vec_indcs <- as.vector(unlist(Map('+', sqnc, list_neigh)))
  #     neighbors[[i]][-vec_indcs] <- NA
  #     neighbors_diss[[i]][-vec_indcs] <- NA
  #   }
  # }

  if (diss_usage == "predictors") {
    if (diss_method %in% c("pca", "pca.nipals", "pls")) {
      if (!is_local) {
        # needs to be scaled/converted on the basis of the scores Xr and Xu
        info_neighbors$diss_xr_xr <- f_diss(euclid_to_mahal(info_neighbors$projection$scores)[1:nrow(Xr), , drop = FALSE],
          diss_method = "euclid",
          center = FALSE, scale = FALSE
        )
      }
      # if (".local" %in% names(neighbor_dots)) {
      #   if (isTRUE(neighbor_dots$local)) {
      #     k_indices <-  apply(info_neighbors$diss_xr_xr,
      #                         MARGIN = 1,
      #                         FUN = function(x, pre_k) order(x)[1:pre_k],
      #                         pre_k = neighbor_dots$pre_k + 1)
      #
      #     info_neighbors$diss_xr_xr <- local_ortho_diss(k_index_matrix = k_indices,
      #                                                   Xr = Xr, Yr = Yr,
      #                                                   Xu = NULL, diss_method = diss_method,
      #                                                   pc_selection = pc_selection,
      #                                                   center = center, scale = scale,
      #                                                   allow_parallel = allow_parallel,
      #                                                   ...)
      #
      #     info_neighbors$local_n_components <- info_neighbors$diss_xr_xr$local_n_components
      #     info_neighbors$diss_xr_xr <- info_neighbors$diss_xr_xr$dissimilarity_m
      #   }
      # }
    } else {
      # needs to be scaled on the basis of the scores Xr and Xu
      info_neighbors$diss_xr_xr <- dissimilarity(scale(rbind(Xr, Xu),
        center = center,
        scale = scale
      )[1:nrow(Xr), ],
      diss_method = diss_method,
      center = FALSE,
      scale = FALSE, ...
      )$dissimilarity
    }
  }

  info_neighbors
}

#' @title Cross validation for PLS regression
#' @description for internal use only!
#' @keywords internal
pls_cv <- function(x, y, ncomp,
                   method = c("pls", "wapls"),
                   center = TRUE, scale,
                   min_component = 1,
                   new_x = matrix(0, 1, 1),
                   weights = NULL,
                   p = 0.75,
                   number = 10,
                   group = NULL,
                   retrieve = TRUE,
                   tune = TRUE,
                   max_iter = 1, tol = 1e-6, 
                   seed = NULL,
                   modified = FALSE) {
  min_allowed <- (floor(min(ncol(x), nrow(x)) * p)) - 1
  
  if (min_allowed < ncomp) {
    ncomp <- min_allowed
  }
  
  if (modified) {
    algorithm <- "mpls"
  } else {
    algorithm <- "pls"
  }
  
  cv_samples <- sample_stratified(
    y = y,
    p = p,
    number = number,
    group = group,
    replacement = FALSE,
    seed = seed
  )
  
  if (method == "wapls" & retrieve & tune) {
    pls_c <- c(min_component, ncomp)
    seq_pls <- min(pls_c):max(pls_c)
    search_grid <- expand.grid(
      min_component = seq_pls,
      max_component = seq_pls,
      KEEP.OUT.ATTRS = FALSE
    )
    search_grid <- as.matrix(search_grid[search_grid$min_component < search_grid$max_component, ])
    rownames(search_grid) <- 1:nrow(search_grid)
    fit_method <- "completewapls1"
  } else {
    search_grid <- matrix(0, 0, 0)
    fit_method <- method
  }
  
  cv_results <- opls_cv_cpp(
    X = x,
    Y = y,
    scale = scale,
    method = fit_method,
    mindices = cv_samples$hold_in,
    pindices = cv_samples$hold_out,
    min_component = min_component,
    ncomp = ncomp,
    new_x = new_x,
    maxiter = max_iter,
    tol = tol,
    wapls_grid = search_grid, 
    algorithm = algorithm
  )
  
  val <- NULL
  val$resamples <- cv_samples$hold_out
  
  if (method == "pls") {
    val$cv_results <- data.table(
      npls = 1:ncomp,
      rmse_cv = rowMeans(cv_results$rmse_seg),
      st_rmse_cv = rowMeans(cv_results$st_rmse_seg),
      rmse_sd_cv = get_col_sds(t(cv_results$rmse_seg)),
      r2_cv = rowMeans(cv_results$rsq_seg)
    )
    best_pls_c <- val$cv_results$npls[which(val$cv_results$rmse_cv == min(val$cv_results$rmse_cv))]
    val$best_pls_c <- best_pls_c
    
    if (retrieve) {
      if (!is.null(weights)) {
        x <- sweep(x, 1, weights, "*")
        y <- y * weights
      }
      if (tune) {
        val$models <- opls_get_basics(
          X = x,
          Y = y,
          scale = scale,
          ncomp = best_pls_c,
          maxiter = max_iter,
          tol = tol, 
          algorithm = algorithm
        )
      } else {
        val$models <- opls_get_basics(
          X = x,
          Y = y,
          scale = scale,
          ncomp = ncomp,
          maxiter = max_iter,
          tol = tol, 
          algorithm = algorithm
        )
      }
    }
  }
  
  if (method == "wapls" & retrieve & !tune) {
    val$compweights <- cv_results$compweights
    val$cv_results <- data.table(
      min_component = min_component,
      max_component = ncomp,
      rmse_cv = mean(cv_results$rmse_seg),
      st_rmse_cv = mean(cv_results$st_rmse_seg),
      rmse_sd_cv = sd(cv_results$rmse_seg),
      r2_cv = mean(cv_results$rsq_seg)
    )
    
    
    if (!is.null(weights)) {
      x <- sweep(x, 1, weights, "*") ### MODIFIED
      y <- y * weights
    }
    val$models <- opls_get_basics(
      X = x,
      Y = y,
      scale = scale,
      ncomp = ncomp,
      maxiter = max_iter,
      tol = tol, 
      algorithm = algorithm
    )
  }
  
  if (method == "wapls" & retrieve & tune) {
    val$cv_results <- data.table(search_grid,
                                 rmse_cv = rowMeans(cv_results$rmse_seg),
                                 st_rmse_cv = rowMeans(cv_results$st_rmse_seg),
                                 rmse_sd_cv = get_col_sds(t(cv_results$rmse_seg)),
                                 r2_cv = rowMeans(cv_results$rsq_seg)
    )
    opmls <- which.min(val$cv_results$rmse_cv)
    val$best_pls_c$min_component <- val$cv_results$min_component[opmls]
    val$best_pls_c$max_component <- val$cv_results$max_component[opmls]
    val$compweights <- cv_results$compweights
    
    
    if (!is.null(weights)) {
      x <- sweep(x, 1, weights, "*")
      y <- y * weights
    }
    
    tmp_compweights <- val$compweights[val$best_pls_c$min_component:val$best_pls_c$max_component]
    val$compweights[] <- 0
    tmp_compweights <- tmp_compweights / sum(tmp_compweights)
    val$compweights[val$best_pls_c$min_component:val$best_pls_c$max_component] <- tmp_compweights
    
    val$models <- opls_get_basics(
      X = x,
      Y = y,
      scale = scale,
      ncomp = ncomp,
      maxiter = max_iter,
      tol = tol, 
      algorithm = algorithm
    )
  }
  val
}


#' @title Internal function for computing the weights of the PLS components
#' necessary for weighted average PLS
#' @description internal
#' @keywords internal
#' @param pls_model either an object returned by the \code{pls_cv} function or an
#' object as returned by the \code{opls_get_basics} function which contains a pls model.
#' @param original_x the original spectral matrix which was used for calibrating the
#' pls model.
#' @param type type of weight to be computed. The only available option (for
#' the moment) is \code{"w1"}. See details on the \code{mbl} function where it
#' is explained how \code{"w1"} is computed whitin the \code{"wapls"}
#' regression.
#' @param new_x a vector of a new spectral observation. When "w1" is selected, new_x
#' must be specified.
#' @param pls_c a vector of length 2 which contains both the minimum and maximum
#' number of PLS components for which the weights must be computed.
#' @return \code{get_wapls_weights} returns a vector of weights for each PLS
#' component specified
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
get_wapls_weights <- function(pls_model, original_x, type = "w1", new_x = NULL, pls_c) {
  if (!is.null(pls_model$models)) {
    rmse_cv <- pls_model$cv_results$rmse_cv
    pls_model <- pls_model$models
    is_plsCv <- TRUE
  } else {
    is_plsCv <- FALSE
  }
  if (max(pls_c) > max(pls_model$ncomp)) {
    stop("the maximum number of PLS components specified in pls_c exceeds the number of PLS components contained in pls_model!")
  }
  if (type == "w1" & is.null(new_x)) {
    stop("new_x is missng! When type = 'w2', a vector of a new spectral sample must be specified")
  }
  if (type == "w1" & (length(new_x) != ncol(original_x))) {
    stop("the length of the vector new_x does not match with the number of variables of original_x")
  }
  if (length(pls_c) != 2) {
    stop("'pls_c' must be a vector of length 2 which specifies both the minimum and the maximum number of PLS components")
  }
  if (pls_c[[1]] > pls_c[[2]]) {
    stop("for the vector of 'pls_c' the first value must be the minimum number of PLS components")
  }

  min_component <- pls_c[[1]]
  max_component <- pls_c[[2]]


  whgt <- get_local_pls_weights(
    projection_mat = pls_model$projection_mat,
    xloadings = pls_model$X_loadings,
    coefficients = pls_model$coefficients,
    new_x = new_x,
    min_component = min_component,
    max_component = max_component,
    scale = ifelse(nrow(pls_model$transf$Xscale) == 1, TRUE, FALSE),
    Xcenter = pls_model$transf$Xcenter,
    Xscale = pls_model$transf$Xscale
  )[1, ]

  return(whgt)
}


#' @title An iterator for local prediction data in mbl
#' @description internal function. It collects only the data necessary to
#' execute a local prediction for the mbl function based on a list of neighbors.
#' Not valid for local dissmilitary (e.g. for ortho_diss(...., .local = TRUE))
#' @param Xr the Xr matrix in mbl.
#' @param Xu the Xu matrix in mbl.
#' @param Yr the Yr matrix in mbl.
#' @param Yu the Yu matrix in mbl.
#' @param diss_usage a character string indicating if the dissimilarity data
#' will be used as predictors ("predictors") or not ("none").
#' @param neighbor_indices a matrix with the indices of neighbors of every Xu
#' found in Xr.
#' @param neighbor_diss a matrix with the dissimilarity socres for the neighbors
#' of every Xu found in Xr. This matrix is organized in the same way as
#' \code{neighbor_indices}.
#' @param diss_xr_xr a dissimilarity matrix between sampes in Xr.
#' @param group a factor representing the group labes of Xr.
#' @return an object of \code{class} iterator giving the following list:
#' itemize{
#' \item{ith_xr:}{ the Xr data of the neighbors for the ith observation (if
#' \code{diss_usage = "predictors"}, this data is combined with the local
#' dissmilarity scores of the neighbors of Xu)}
#' \item{ith_yr:}{ the Yr data of the neighbors for the ith observation}
#' \item{ith_xu:}{ the ith Xu observation (if \code{diss_usage = "predictors"},
#' this data is combined with the local dissmilarity scores to its Xr neighbors}
#' \item{ith_yu:}{ the ith Yu observation}
#' \item{ith_neigh_diss:}{ the dissimilarity scores of the neighbors for the ith
#' observation}
#' \item{ith_group:}{ the group labels for ith_xr}
#' \item{n_k:}{ the number of neighbors}
#' }
#' @details isubset will look at the order of knn in each col of D and
#' re-organize the rows of x accordingly
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
ith_mbl_neighbor <- function(Xr, Xu, Yr, Yu = NULL,
                             diss_usage = "none",
                             neighbor_indices,
                             neighbor_diss = NULL,
                             diss_xr_xr = NULL,
                             group = NULL) {
  k_neighbors <- colSums(!is.na(neighbor_indices))
  iter_k_neighbors <- iter(k_neighbors, by = "cell")
  iter_neighbors <- iter(neighbor_indices, by = "col")
  iter_xu <- iter(Xu, by = "row")
  iter_yu <- iter(Yu, by = "cell")
  neighbor_diss <- t(neighbor_diss)
  iter_xr_xu_diss <- iter(neighbor_diss, by = "row")
  group

  nextEl <- function() {
    ith_neighborhood <- nextElem(iter_neighbors)
    ith_ks <- 1:nextElem(iter_k_neighbors)
    ith_neighborhood <- ith_neighborhood[ith_ks]

    ith_xr_neighbors <- Xr[ith_neighborhood, , drop = FALSE]

    ith_yr_neighbors <- Yr[ith_neighborhood, , drop = FALSE]
    ith_xu <- nextElem(iter_xu)
    ith_neigh_diss <- nextElem(iter_xr_xu_diss)[, ith_ks, drop = FALSE]

    if (!is.null(Yu)) {
      ith_yu <- nextElem(iter_yu)
    } else {
      ith_yu <- NULL
    }

    if (diss_usage == "predictors") {
      ith_local_xr_xr_diss <- diss_xr_xr[ith_neighborhood, ith_neighborhood]
      colnames(ith_local_xr_xr_diss) <- paste0("k_", ith_ks)
      ith_xr_neighbors <- cbind(ith_local_xr_xr_diss, ith_xr_neighbors)
      ith_xu <- cbind(ith_neigh_diss, ith_xu)
    }

    if (!is.null(group)) {
      ith_group <- factor(group[ith_neighborhood])
    } else {
      ith_group <- NULL
    }

    it_neigh <- list(
      ith_xr = ith_xr_neighbors,
      ith_yr = ith_yr_neighbors,
      ith_xu = ith_xu,
      ith_yu = ith_yu,
      ith_neig_indices = ith_neighborhood,
      ith_neigh_diss = ith_neigh_diss,
      local_index_nearest = which.min(ith_neigh_diss)[1],
      ith_group = ith_group,
      n_k = max(ith_ks)
    )
    it_neigh
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}




#' @title A function to obtain the local neighbors based on dissimilarity
#' matrices from orthogonal projections.
#' @description internal function. This function is used to obtain the local
#' neighbors based on dissimilarity matrices from orthogonal projections. These
#' neighbors are obatin from an orthogonal projection on a set of precomputed
#' neighbors. This function is used internally by the mbl fucntion.
#' ortho_diss(, .local = TRUE) operates in the same way, however for mbl, it is
#' more efficient to do the re-search of the neighbors inside its main for loop
#'
#' @param ith_xr the set of neighbors of a Xu observation found in Xr
#' @param ith_xu the Xu observation
#' @param ith_yr the response values of the set of neighbors of the Xu
#' observation found in Xr
#' @param ith_yu the response value of the xu observation
#' @param diss_usage a character string indicating if the dissimilarity data
#' will be used as predictors ("predictors") or not ("none").
#' @param ith_neig_indices a vector of the original indices of the Xr neighbors.
#' @param k the number of nearest neighbors to select from the already
#' identified neighbors
#' @param k_diss the distance threshold to select the neighbors from the already
#' identified neighbors
#' @param k_range a min and max  number of allowed neighbors when \code{k_diss}
#' is used
#' @param spike a vector with the indices of the observations forced to be
#' retained as neighbors. They have to be present in all the neighborhoods and
#' at the top of \code{neighbor_indices}.
#' @param diss_method the ortho_diss() method
#' @param pc_selection the pc_selection argument as in ortho_diss()
#' @param ith_group the vector containing the group labes of \code{ith_xr}.
#' @param center center the data in the local diss computation?
#' @param scale scale the data in the local diss computation?
#' @return a list:
#' itemize{
#' \item{ith_xr:}{ the new Xr data of the neighbors for the ith observation (if
#' \code{diss_usage = "predictors"}, this data is combined with the local
#' dissmilarity scores of the neighbors of Xu)}
#' \item{ith_yr:}{ the new Yr data of the neighbors for the ith observation}
#' \item{ith_xu:}{ the ith Xu observation (if \code{diss_usage = "predictors"},
#' this data is combined with the local dissmilarity scores to its Xr neighbors}
#' \item{ith_yu:}{ the ith Yu observation}
#' \item{ith_neigh_diss:}{ the new dissimilarity scores of the neighbors for the ith
#' observation}
#' \item{ith_group:}{ the group labels for the new ith_xr}
#' \item{n_k:}{ the number of neighbors}
#' \item{ith_components:}{ the number of components used}
#' }
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#'
get_ith_local_neighbors <- function(ith_xr, ith_xu, ith_yr, ith_yu = NULL,
                                    diss_usage = "none",
                                    ith_neig_indices,
                                    k = NULL, k_diss = NULL, k_range = NULL,
                                    spike = NULL, diss_method, pc_selection,
                                    ith_group = NULL,
                                    center, scale, ...) {
  is_predictors <- diss_usage == "predictors"

  neighbor_dots <- list(...)
  # local dissimilarities option is always
  # deactivated, as this will be calculated inside the local loops for
  # modeling. This Saves memory
  is_local <- isTRUE(neighbor_dots$.local)
  if (is_local) {
    neighbor_dots$.local <- FALSE
  }

  if (!is.null(spike)) {
    # 1:length(spike) is used in spike as the neighbors are reorganized and
    # and moved to the top of the Xr data. So the first length(spike) are
    # retained
    reorg_spike <- 1:length(spike)
  } else {
    reorg_spike <- NULL
  }

  # use do.call to be able to change loal to FALSE in ... (if provided)
  local_diss <- do.call(
    ortho_diss,
    c(
      list(
        Xr = ith_xr, Xu = ith_xu,
        Yr = ith_yr,
        pc_selection = pc_selection,
        diss_method = diss_method,
        center = center,
        scale = scale,
        compute_all = is_predictors,
        return_projection = FALSE,
        allow_parallel = FALSE
      ),
      neighbor_dots
    )
  )



  if (is_predictors) {
    indx_xu <- ncol(local_diss$dissimilarity)

    loc_neigh <- diss_to_neighbors(
      local_diss$dissimilarity[-indx_xu, indx_xu, drop = FALSE],
      k = k, k_diss = k_diss, k_range = k_range,
      spike = reorg_spike,
      return_dissimilarity = FALSE
    )

    neigh_indcs <- loc_neigh$neighbors
    ith_local_xr_xr_diss <- local_diss$dissimilarity[neigh_indcs, neigh_indcs]
    colnames(ith_local_xr_xr_diss) <- paste0(
      "k_",
      1:ncol(ith_local_xr_xr_diss)
    )
    ith_xr <- cbind(
      ith_local_xr_xr_diss,
      ith_xr[loc_neigh$neighbors, ]
    )
    ith_yr <- ith_yr[neigh_indcs, , drop = FALSE]
    ith_xu <- cbind(t(loc_neigh$neighbors_diss), ith_xu)
  } else {
    loc_neigh <- diss_to_neighbors(local_diss$dissimilarity,
      k = k, k_diss = k_diss, k_range = k_range,
      spike = reorg_spike,
      return_dissimilarity = FALSE
    )

    ith_xr <- ith_xr[loc_neigh$neighbors, ]
    ith_yr <- ith_yr[loc_neigh$neighbors, , drop = FALSE]
  }
  loc_neighbors <- ith_neig_indices[loc_neigh$neighbors]
  ith_neigh_diss <- loc_neigh$neighbors_diss

  if (!is.null(ith_group)) {
    ith_group <- ith_group[loc_neigh$neighbors]
  }


  it_neigh <- list(
    ith_xr = ith_xr,
    ith_yr = ith_yr,
    ith_xu = ith_xu,
    ith_yu = ith_yu,
    ith_neig_indices = loc_neighbors,
    ith_neigh_diss = ith_neigh_diss,
    local_index_nearest = which.min(ith_neigh_diss)[1],
    ith_group = ith_group,
    n_k = length(loc_neighbors),
    ith_components = local_diss$n_components
  )
  it_neigh
}

#' @title Standard deviation of columns
#' @description For internal use only!
#' @keywords internal
get_col_sds <- function(x) {
  return(as.vector(get_column_sds(x)))
}


#' @title Local multivariate regression
#' @description internal
#' @keywords internal
fit_and_predict <- function(x, y, pred_method, scale = FALSE, weights = NULL,
                            newdata, pls_c = NULL, CV = FALSE,
                            tune = FALSE, number = 10, p = 0.75,
                            group = NULL, noise_variance = 0.001,
                            range_prediction_limits = TRUE,
                            pls_max_iter = 1, pls_tol = 1e-6, modified = FALSE, seed = NULL) {
  if (is.null(weights)) {
    weights <- 1
  }

  if (any(get_col_sds(x) == 0)) {
    warning("One of the variables has zero variance. Data will not be scale")
    scale <- FALSE
  }

  if (modified) {
    algorithm <- "mpls"
  } else {
    algorithm <- "pls"
  }
    

  results <- cv_val <- pred <- NULL
  if (pred_method == "gpr") {
    if (CV) {
      cv_val <- gaussian_pr_cv(
        x = x,
        y = y,
        scale = scale,
        weights = weights,
        p = p,
        number = number,
        group = group,
        retrieve = "final_model", 
        seed = seed
      )
      fit <- cv_val$model
    } else {
      x <- sweep(x, 1, weights, "*")
      y <- y * weights
      fit <- gaussian_process(
        X = x,
        Y = y,
        noisev = noise_variance,
        scale = scale
      )
    }

    pred <- predict_gaussian_process(
      Xz = fit$Xz,
      alpha = fit$alpha,
      newdata = newdata,
      scale = fit$is_scaled,
      Xcenter = fit$Xcenter,
      Xscale = fit$Xscale,
      Ycenter = fit$Ycenter,
      Yscale = fit$Yscale
    )[[1]]
  }
  if (pred_method == "pls") {
    if (CV) {
      cv_val <- pls_cv(
        x = x,
        y = y,
        ncomp = pls_c,
        method = "pls",
        center = TRUE,
        scale = scale,
        weights = weights,
        p = p,
        number = number,
        group = group,
        retrieve = TRUE,
        tune = tune,
        max_iter = pls_max_iter,
        tol = pls_tol, 
        modified = modified, 
        seed = seed
      )
      fit <- cv_val$models
      ncomp <- cv_val$best_pls_c
    } else {
      x <- sweep(x, 1, weights, "*")
      y <- y * weights
      fit <- opls_get_basics(
        X = x,
        Y = y,
        scale = scale,
        ncomp = pls_c,
        maxiter = pls_max_iter,
        tol = pls_tol, 
        algorithm = algorithm
      )
      ncomp <- pls_c
    }
    pred <- predict_opls(
      bo = fit$bo,
      b = fit$coefficients,
      ncomp = ncomp,
      newdata = newdata,
      scale = scale,
      Xscale = fit$transf$Xscale
    )[, ncomp]
  }

  if (pred_method == "wapls") {
    if (CV) {
      cv_val <- pls_cv(
        x = x, y = y,
        ncomp = pls_c[[2]],
        method = "wapls",
        center = TRUE,
        scale = scale,
        min_component = pls_c[[1]],
        new_x = newdata,
        weights = weights,
        p = p,
        number = number,
        group = group,
        retrieve = TRUE,
        tune = tune,
        max_iter = pls_max_iter,
        tol = pls_tol, 
        modified = modified, 
        seed = seed
      )

      fit <- cv_val$models
      if (!tune) {
        cv_val$cv_results <- data.table(
          min_component = pls_c[[1]],
          max_component = pls_c[[2]],
          rmse_cv = cv_val$cv_results$rmse_cv,
          st_rmse_cv = cv_val$cv_results$st_rmse_cv,
          rmse_sd_cv = cv_val$cv_results$rmse_sd_cv,
          r2_cv = cv_val$cv_results$r2_cv
        )
      }


      ## here all the weights are output from 1 to plsMax
      w <- cv_val$compweights
    } else {
      x <- sweep(x, 1, weights, "*")
      y <- y * weights
      fit <- opls_get_basics(
        X = x,
        Y = y,
        scale = scale,
        ncomp = pls_c[[2]],
        maxiter = pls_max_iter,
        tol = pls_tol, 
        algorithm = algorithm
      )

      # compute weights for PLS components selected (from plsMin to plsMax)
      w <- rep(0, pls_c[[2]])
      w[pls_c[[1]]:pls_c[[2]]] <- get_wapls_weights(
        pls_model = fit,
        original_x = x,
        type = "w1",
        new_x = newdata,
        pls_c = pls_c
      )
    }

    # compute the weighted average of the multiple PLS predictions
    new_x_prediction <- predict_opls(
      bo = fit$bo,
      b = fit$coefficients,
      ncomp = pls_c[[2]],
      newdata = newdata,
      scale = scale,
      Xscale = fit$transf$Xscale
    )
    pred <- sum(new_x_prediction * w) # weighted preds
  }

  if (range_prediction_limits) {
    if (pred > max(y)) {
      pred <- max(y)
    }
    if (pred < min(y)) {
      pred <- min(y)
    }
  }

  results$validation <- cv_val
  results$prediction <- pred

  results
}

#' @title Cross validation for Gaussian process regression
#' @description internal
#' @keywords internal
gaussian_pr_cv <- function(x,
                           y,
                           scale,
                           weights = NULL,
                           p = 0.75,
                           number = 10,
                           group = NULL,
                           noise_variance = 0.001,
                           retrieve = c("final_model", "none"), 
                           seed = NULL) {

  ## Create the resampling groups
  cv_samples <- sample_stratified(
    y = y,
    p = p,
    number = number,
    group = group,
    replacement = FALSE, 
    seed = seed
  )

  cv_results <- gaussian_process_cv(
    X = x,
    Y = y,
    mindices = cv_samples$hold_in,
    pindices = cv_samples$hold_out,
    noisev = noise_variance,
    scale = scale
  )

  val <- NULL
  val$resamples <- cv_samples$hold_out
  val$cv_results <- data.table(
    rmse_cv = mean(cv_results$rmse_seg),
    st_rmse_cv = mean(cv_results$st_rmse_seg),
    rmse_sd_cv = sd(cv_results$rmse_seg),
    r2_cv = mean(cv_results$rsq_seg)
  )

  if (retrieve == "final_model") {
    if (is.null(weights)) {
      x <- sweep(x, 1, weights, "*")
      y <- y * weights
    }
    val$model <- gaussian_process(
      X = x,
      Y = y,
      noisev = noise_variance,
      scale = scale
    )
  }
  val
}
