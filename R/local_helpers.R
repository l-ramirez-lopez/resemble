#' @title Get neighbor information
#' @description Gathers neighborhood information for Xu observations in Xr.
#' Required during local regressions.
#' @details
#' For local PCA and PLS distances, local dissimilarity matrices are not
#' computed here (cheaper to compute during local regressions). Instead,
#' global distances are output for later local dissimilarity computation.
#' @keywords internal
#' @noRd
get_neighbor_info <- function(
    Xr, Xu, 
    diss_method,
    Yr = NULL,
    neighbors,
    spike = NULL,
    return_dissimilarity = FALSE,
    diss_usage = "none"
) {
  
  n_xr <- nrow(Xr)
  
  # Determine if using local dissimilarities
  is_local <- inherits(diss_method, "diss_local_pca") || 
    inherits(diss_method, "diss_local_pls")
  
  if (is_local) {
    # For local methods, use pre_k for initial neighbor search
    neighbors <- neighbors_k(diss_method$pre_k)
    return_diss <- FALSE
  } else {
    return_diss <- return_dissimilarity
  }
  
  # Search neighbors
  info_neighbors <- search_neighbors(
    Xr = Xr, 
    Xu = Xu, 
    diss_method = diss_method,
    Yr = Yr,
    neighbors = neighbors,
    spike = spike,
    return_dissimilarity = return_diss
  )
  
  # Compute Xr-Xr dissimilarity if needed for predictors
  if (diss_usage == "predictors" && !is_local) {
    is_ortho <- inherits(diss_method, "diss_pca") || 
      inherits(diss_method, "diss_pls")
    if (is_ortho) {
      # Convert scores to Mahalanobis space, compute Euclidean on Xr portion
      scores_mahal <- euclid_to_mahal(info_neighbors$projection$scores)
      info_neighbors$diss_xr_xr <- f_diss(
        scores_mahal[1:n_xr, , drop = FALSE],
        diss_method = "euclid",
        center = FALSE, 
        scale = FALSE
      )
    } else {
      # Non-orthogonal methods: compute directly on scaled data
      center <- diss_method$center %||% TRUE
      scale <- diss_method$scale %||% FALSE
      
      Xr_scaled <- scale(
        rbind(Xr, Xu),
        center = center,
        scale = scale
      )[1:n_xr, , drop = FALSE]
      
      info_neighbors$diss_xr_xr <- dissimilarity(
        Xr = Xr_scaled,
        diss_method = diss_method
      )$dissimilarity
    }
  }
  
  info_neighbors
}


#' @title Cross validation for PLS regression
#' @description for internal use only!
#' @keywords internal
pls_cv <- function(
    x, y, ncomp,
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
    algorithm = c("pls", "mpls", "simpls")
) {
  min_allowed <- (floor(min(ncol(x), nrow(x)) * p)) - 1
  
  if (min_allowed < ncomp) {
    ncomp <- min_allowed
  }
  
  modified <- ifelse(algorithm == "mpls", TRUE, FALSE)

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
    val$cv_results <- data.frame(
      ncomp = 1:ncomp,
      rmse_cv = rowMeans(cv_results$rmse_seg),
      st_rmse_cv = rowMeans(cv_results$st_rmse_seg),
      rmse_sd_cv = get_col_sds(t(cv_results$rmse_seg)),
      r2_cv = rowMeans(cv_results$rsq_seg)
    )
    best_pls_c <- val$cv_results$ncomp[which(val$cv_results$rmse_cv == min(val$cv_results$rmse_cv))]
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
    val$cv_results <- data.frame(
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
    val$cv_results <- data.frame(search_grid,
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
#' is explained how \code{"w1"} is computed within the \code{"wapls"}
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
  
  do_scaling <- ifelse(nrow(pls_model$transf$Xscale) == 1, TRUE, FALSE)
  do_centering <- ifelse(nrow(pls_model$transf$Xcenter) == 1, TRUE, FALSE)
  
  if (do_scaling) {
    scaling_vec <- pls_model$transf$Xscale[1, ]
  } else {
    scaling_vec <- matrix(NA, 0, 0)
  }
  
  if (do_centering) {
    centering_vec <- pls_model$transf$Xcenter[1, ]
  } else {
    scaling_vec <- matrix(NA, 0, 0)
  }
  
  
  whgt <- get_local_pls_weights(
    projection_mat = pls_model$projection_mat,
    xloadings = pls_model$X_loadings,
    coefficients = pls_model$coefficients,
    new_x = new_x,
    min_component = min_component,
    max_component = max_component,
    scale = do_scaling,
    Xcenter = centering_vec,
    Xscale = scaling_vec
  )[1, ]
  
  return(whgt)
}


#' @title An iterator for local prediction data in mbl
#' @description internal function. It collects only the data necessary to
#' execute a local prediction for the mbl function based on a list of neighbors.
#' Not valid for local dissmilitary (e.g. for ortho_diss(...., .local = TRUE))
#' @param Xr the Xr matrix in mbl.
#' @param Xu the Xu matrix in mbl. Default \code{NULL}. If not provided, the 
#' function will iterate for each \code{\{Yr, Xr\}} to get the respective neighbors.
#' @param Yr the Yr matrix in mbl.
#' @param Yu the Yu matrix in mbl. Default \code{NULL}. 
#' @param diss_usage a character string indicating if the dissimilarity data
#' will be used as predictors ("predictors") or not ("none").
#' @param neighbor_indices a matrix with the indices of neighbors of every Xu
#' found in Xr.
#' @param neighbor_diss a matrix with the dissimilarity socres for the neighbors
#' of every Xu found in Xr. This matrix is organized in the same way as
#' \code{neighbor_indices}.
#' @param diss_xr_xr a dissimilarity matrix between sampes in Xr.
#' @param group a factor representing the group labels of Xr.
#' @return an object of \code{class} iterator giving the following list:
#' \itemize{
#' \item{ith_xr: the Xr data of the neighbors for the ith observation (if
#' \code{diss_usage = "predictors"}, this data is combined with the local
#' dissmilarity scores of the neighbors of Xu (or Xr if Xu was not provided))}
#' \item{ith_yr: the Yr data of the neighbors for the ith observation}
#' \item{ith_xu: the ith Xu observation (or Xr if Xu was not provided). 
#' If \code{diss_usage = "predictors"}, this data is combined with the local 
#' dissmilarity scores to its Xr neighbors.}
#' \item{ith_yu: the ith Yu observation (or Yr observation if Xu was not provided).}
#' \item{ith_neigh_diss: the dissimilarity scores of the neighbors for the ith
#' observation.}
#' \item{ith_group: the group labels for ith_xr.}
#' \item{n_k: the number of neighbors.}
#' }
#' @details isubset will look at the order of knn in each col of D and
#' re-organize the rows of x accordingly
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
ith_mbl_neighbor <- function(
    Xr, Xu = NULL, Yr, Yu = NULL,
    diss_usage = "none",
    neighbor_indices,
    neighbor_diss = NULL,
    diss_xr_xr = NULL,
    group = NULL
) {
  k_neighbors <- colSums(!is.na(neighbor_indices))
  iter_k_neighbors <- iter(k_neighbors, by = "cell", recycle=T)
  iter_neighbors <- iter(neighbor_indices, by = "col", recycle=T)
  
  if (is.null(Xu)) {
    iter_xu <- iter(Xr, by = "row" , recycle=T)
    iter_yu <- iter(Yr, by = "cell", recycle=T)
    iterate_yu <- TRUE
  } else {
    iter_xu <- iter(Xu, by = "row", recycle=T)
    iter_yu <- iter(Yu, by = "cell", recycle=T)
    iterate_yu <- !is.null(Yu)
  }
  neighbor_diss <- t(neighbor_diss)
  iter_xr_xu_diss <- iter(neighbor_diss, by = "row", recycle=T)
  group
  nextEl <- function() {
    ith_neighborhood <- nextElem(iter_neighbors)
    ith_ks <- 1:nextElem(iter_k_neighbors)
    ith_neighborhood <- ith_neighborhood[ith_ks]
    
    ith_xr_neighbors <- Xr[ith_neighborhood, , drop = FALSE]
    
    ith_yr_neighbors <- Yr[ith_neighborhood, , drop = FALSE]
    ith_xu <- nextElem(iter_xu)
    ith_neigh_diss <- nextElem(iter_xr_xu_diss)[, ith_ks, drop = FALSE]
    
    if (iterate_yu) {
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


#' @title Get local neighbors based on orthogonal dissimilarity
#' @description Internal function for obtaining local neighbors based on 
#' dissimilarity matrices from orthogonal projections. These neighbors are 
#' obtained from an orthogonal projection on a set of precomputed neighbors.
#' Used internally by \code{mbl()}.
#' @param ith_xr The set of neighbors of a Xu observation found in Xr.
#' @param ith_xu The Xu observation.
#' @param ith_yr The response values of the neighbors.
#' @param ith_yu The response value of the Xu observation.
#' @param diss_usage Character indicating if dissimilarity data will be used 
#'   as predictors ("predictors") or not ("none").
#' @param ith_neig_indices Vector of original indices of the Xr neighbors.
#' @param neighbors A neighbor selection object from \code{neighbors_k()} or 
#'   \code{neighbors_diss()}.
#' @param spike Vector of observation indices forced to be retained as neighbors.
#' @param diss_method A local dissimilarity method object (\code{diss_local_pca()} 
#'   or \code{diss_local_pls()}).
#' @param ith_group Vector containing group labels for \code{ith_xr}.
#' @param mbl_is_parallel Logical indicating if \code{mbl()} is running in 
#'   parallel. If \code{TRUE}, nested parallelism is prevented by forcing 
#'   \code{allow_parallel = FALSE} in the dissimilarity computation.
#' @return A list containing:
#' \itemize{
#'   \item{ith_xr: Xr data of neighbors (combined with local dissimilarity 
#'     scores if \code{diss_usage = "predictors"})}
#'   \item{ith_yr: Yr data of neighbors}
#'   \item{ith_xu: Xu observation (combined with local dissimilarity scores 
#'     if \code{diss_usage = "predictors"})}
#'   \item{ith_yu: Yu observation}
#'   \item{ith_neig_indices: Original indices of selected neighbors}
#'   \item{ith_neigh_diss: Dissimilarity scores of neighbors}
#'   \item{local_index_nearest: Index of nearest neighbor}
#'   \item{ith_group: Group labels for ith_xr}
#'   \item{n_k: Number of neighbors}
#'   \item{ith_ncomp: Number of components used}
#' }
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @keywords internal
#' @noRd
get_ith_local_neighbors <- function(
    ith_xr, ith_xu, ith_yr, ith_yu = NULL,
    diss_usage = "none",
    ith_neig_indices,
    neighbors,
    spike = NULL, 
    diss_method,
    ith_group = NULL,
    mbl_is_parallel = FALSE
) {
  is_predictors <- diss_usage == "predictors"
  
  # Prevent nested parallelism if mbl is running in parallel
  if (mbl_is_parallel) {
    diss_method$allow_parallel <- FALSE
  }
  
  # Handle spike indices (only positive values are forced in)
  if (!is.null(spike)) {
    n_spike_hold <- sum(spike > 0)
    if (n_spike_hold > 0) {
      reorg_spike <- seq_len(n_spike_hold)
    } else {
      reorg_spike <- NULL
    }
  } else {
    reorg_spike <- NULL
  }
  
  # Convert neighbors to internal format
  if (inherits(neighbors, "neighbors_k")) {
    k <- max(neighbors$k)
    k_diss <- NULL
    k_range <- NULL
  } else {
    k <- NULL
    k_diss <- max(neighbors$threshold)
    k_range <- c(neighbors$k_min, neighbors$k_max)
  }
  
  # Compute local dissimilarity
  local_diss <- ortho_diss(
    Xr = ith_xr, 
    Xu = ith_xu,
    Yr = ith_yr,
    diss_method = diss_method,
    compute_all = is_predictors,
    return_projection = FALSE
  )
  
  if (is_predictors) {
    indx_xu <- ncol(local_diss$dissimilarity)
    
    loc_neigh <- diss_to_neighbors(
      local_diss$dissimilarity[-indx_xu, indx_xu, drop = FALSE],
      k = k, 
      k_diss = k_diss, 
      k_range = k_range,
      spike = reorg_spike,
      return_dissimilarity = FALSE
    )
    
    neigh_indcs <- loc_neigh$neighbors
    ith_local_xr_xr_diss <- local_diss$dissimilarity[neigh_indcs, neigh_indcs]
    colnames(ith_local_xr_xr_diss) <- paste0("k_", seq_len(ncol(ith_local_xr_xr_diss)))
    
    ith_xr <- cbind(ith_local_xr_xr_diss, ith_xr[neigh_indcs, ])
    ith_yr <- ith_yr[neigh_indcs, , drop = FALSE]
    ith_xu <- cbind(t(loc_neigh$neighbors_diss), ith_xu)
  } else {
    loc_neigh <- diss_to_neighbors(
      local_diss$dissimilarity,
      k = k, 
      k_diss = k_diss, 
      k_range = k_range,
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
  
  list(
    ith_xr = ith_xr,
    ith_yr = ith_yr,
    ith_xu = ith_xu,
    ith_yu = ith_yu,
    ith_neig_indices = loc_neighbors,
    ith_neigh_diss = ith_neigh_diss,
    local_index_nearest = which.min(ith_neigh_diss)[1],
    ith_group = ith_group,
    n_k = length(loc_neighbors),
    ith_ncomp = local_diss$ncomp
  )
}


#' @title Standard deviation of columns
#' @description For internal use only!
#' @keywords internal
#' @noRd
get_col_sds <- function(x) {
  return(as.vector(get_column_sds(x)))
}


#' @title Local multivariate regression
#' @description internal
#' @keywords internal
#' @noRd
fit_and_predict <- function(
    x, y, pred_method, scale = FALSE, weights = NULL,
    newdata, pls_c = NULL, CV = FALSE,
    tune = FALSE, number = 10, p = 0.75,
    group = NULL, noise_variance = 0.001,
    range_prediction_limits = TRUE,
    pls_max_iter = 1, pls_tol = 1e-6, algorithm, seed = NULL
) {
  if (is.null(weights)) {
    weights <- 1
  }
  
  if (any(get_col_sds(x) == 0)) {
    warning("Variables with zero variance. Data will not be scaled")
    scale <- FALSE
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
        algorithm = algorithm,
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
        algorithm = algorithm,
        seed = seed
      )
      
      fit <- cv_val$models
      if (!tune) {
        cv_val$cv_results <- data.frame(
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
gaussian_pr_cv <- function(
    x,
    y,
    scale,
    weights = NULL,
    p = 0.75,
    number = 10,
    group = NULL,
    noise_variance = 0.001,
    retrieve = c("final_model", "none"),
    seed = NULL
) {
  
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
  val$cv_results <- data.frame(
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


# =============================================================================
# GH distance computation
# =============================================================================

.compute_gh <- function(proj, n_xr) {
  scores <- proj$scores
  
  # Compute centroid from Xr scores only
  scores_center <- colMeans(scores[seq_len(n_xr), , drop = FALSE])
  
  # Mahalanobis distance from center
  gh_diss <- f_diss(
    Xr = scores,
    Xu = matrix(scores_center, nrow = 1),
    diss_method = "mahalanobis",
    center = FALSE,
    scale = FALSE
  )
  gh_vec <- as.vector(gh_diss)
  
  n_total <- nrow(scores)
  
  result <- list(gh_Xr = gh_vec[seq_len(n_xr)])
  
  if (n_total > n_xr) {
    result$gh_Xu <- gh_vec[(n_xr + 1L):n_total]
  }
  
  result
}

#' @title Cross-validated PLS fitting for gesearch
#' 
#' @description
#' Internal function that fits a PLS model with cross-validation for use in
#' `gesearch()`. Handles both `fit_pls()` and `fit_wapls()` methods.
#'
#' @param X Numeric matrix of predictor variables (observations in rows).
#' @param Y Numeric matrix or vector of response values.
#' @param indices Optional integer vector specifying row indices to use for
#'   training. If `NULL`, all rows are used.
#' @param fit_method A `fit_method` object created by [fit_pls()] or
#'   [fit_wapls()] specifying the PLS fitting parameters.
#' @param number Integer. Number of cross-validation iterations.
#' @param group Optional factor assigning group labels to observations for
#'   leave-group-out cross-validation.
#' @param retrieve Logical. If `TRUE`, return resampling indices.
#' @param tune Logical. If `TRUE`, tune the number of PLS components via
#'   cross-validation.
#' @param seed Optional integer for reproducible resampling.
#' @param p Numeric. Proportion of observations to retain in each
#'   cross-validation fold.
#'
#' @return An object of class `"plsr"` containing:
#' \describe{
#'   \item{models}{Fitted PLS model from `opls_get_basics()`.}
#'   \item{cv_results}{Cross-validation results.}
#'   \item{train_resamples}{Resampling indices (if `retrieve = TRUE`).}
#'   \item{fit_method}{The `fit_method` object used.}
#'   \item{Xscale}{Scaling parameters for predictors.}
#'   \item{call}{The matched call.}
#' }
#'
#' @seealso [gesearch()], [fit_pls()], [fit_wapls()]
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
get_plsr <- function(
    X, Y,
    indices = NULL,
    fit_method,
    number = NULL,
    group = NULL,
    retrieve = FALSE,
    tune = TRUE,
    seed = NULL,
    p = NULL
) {
  
  if (!inherits(fit_method, "fit_method")) {
    stop("'fit_method' must be a fit_*() object", call. = FALSE)
  }
  
  # Ensure matrices
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Subset rows for training
  if (is.null(indices)) {
    idx <- seq_len(nrow(X))
  } else {
    idx <- indices
  }
  
  # Align optional grouping to the selected rows
  grp <- if (!is.null(group)) group[idx] else NULL
  
  # Determine ncomp range based on fit_method type
  if (inherits(fit_method, "fit_wapls")) {
    min_ncomp <- fit_method$min_ncomp
    max_ncomp <- fit_method$max_ncomp
  } else if (inherits(fit_method, "fit_pls")) {
    min_ncomp <- 1L
    max_ncomp <- fit_method$ncomp
  } else {
    stop("'fit_method' must be fit_pls() or fit_wapls()", call. = FALSE)
  }
  
  # ---- Cross-validated PLS -------------------------------------------------
  pls_args <- list(
    x = X[idx, , drop = FALSE],
    y = Y[idx, , drop = FALSE],
    ncomp = max_ncomp,
    method = fit_method$fit_method,
    center = fit_method$center,
    scale = fit_method$scale,
    min_component = min_ncomp,
    number = number,
    group = grp,
    retrieve = retrieve,
    tune = tune,
    max_iter = fit_method$max_iter,
    tol = fit_method$tol,
    seed = seed,
    algorithm = fit_method$method
  )
  
  if (!is.null(p)) pls_args$p <- p
  
  plsval <- do.call(pls_cv, pls_args)
  
  # ---- Fit model for prediction --------------------------------------------
  plsval$models <- opls_get_basics(
    X = X[idx, , drop = FALSE],
    Y = Y[idx, , drop = FALSE],
    ncomp = max_ncomp,
    scale = fit_method$scale,
    maxiter = fit_method$max_iter,
    tol = fit_method$tol, 
    algorithm = fit_method$method
  )
  
  # ---- Map resample indices back to original row indices -------------------
  if (!is.null(indices) && !is.null(plsval$resamples)) {
    rnms <- rownames(plsval$resamples)
    plsval$resamples <- apply(plsval$resamples, 2L, function(r) indices[r])
    rownames(plsval$resamples) <- rnms
  }
 
  # ---- Assemble output -----------------------------------------------------
  obj <- list(
    models = plsval$models,
    cv_results = plsval$cv_results,
    train_resamples = plsval$resamples,
    fit_method = fit_method,
    Xscale = plsval$models$transf$Xscale,
    call = match.call()
  )
  class(obj) <- "plsr"
  obj
}



#' @title Predict from a plsr model (internal)
#' 
#' @description
#' Internal S3 method. Generates predictions for new data using a `"plsr"`
#' object (coefficients and intercepts produced by `opls_get_basics()`).
#'
#' @param object A `"plsr"` object returned by [get_plsr()].
#' @param newdata Numeric matrix or data frame of predictors.
#' @param ncomp Integer specifying the number of components to use. Defaults
#'   to the maximum number available in the model.
#' @param ... Unused. Reserved for future extensions.
#'
#' @return A numeric matrix of predictions with columns named
#'   `"pls1"`, `"pls2"`, etc.
#'
#' @details
#' Delegates to `predict_opls()`. Row names of `newdata` are preserved in
#' the output if present.
#'
#' @seealso [get_plsr()], [fit_pls()], [fit_wapls()]
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#' @noRd
predict.plsr <- function(object, newdata, ncomp = NULL, ...) {
  if (!inherits(object, "plsr")) {
    stop("'object' must be of class 'plsr'", call. = FALSE)
  }
  
  Xnew <- as.matrix(newdata)
  
  # Determine ncomp from fit_method if not provided
  if (is.null(ncomp)) {
    if (inherits(object$fit_method, "fit_wapls")) {
      ncomp <- object$fit_method$max_ncomp
    } else if (inherits(object$fit_method, "fit_pls")) {
      ncomp <- object$fit_method$ncomp
    } else {
      ncomp <- ncol(object$models$coefficients)
    }
  }
  
  yhat <- predict_opls(
    bo = object$models$bo,
    b = object$models$coefficients,
    newdata = Xnew,
    ncomp = ncomp,
    scale = object$fit_method$scale,
    Xscale = object$Xscale
  )
  
  colnames(yhat) <- paste0("ncomp", seq_len(ncol(yhat)))
  rownames(yhat) <- if (!is.null(rownames(newdata))) {
    rownames(newdata)
  } else {
    seq_len(nrow(yhat))
  }
  
  yhat
}

#' @title Validation metrics for PLS predictions (internal)
#' @description Internal helper to compute per-component metrics
#'   (RMSE, R², mean error) comparing predicted columns vs a single
#'   observed response.
#'
#' @param predicted Matrix/data frame of predictions
#'   (\code{n} rows × \code{k} components).
#' @param observed Numeric vector or single-column matrix with observed values.
#'
#' @return A data frame with columns:
#'   \itemize{
#'   \item \code{npls}: component index (1..\eqn{k})
#'   \item \code{r2}: squared correlation between predictions and observed
#'   \item \code{rmse}: root mean squared error
#'   \item \code{me}: mean error (signed bias)
#'   }
#'
#' @details Only complete cases in \code{observed} are used. The residuals are
#'   computed as \eqn{\hat{y}_{i,k} - y_i}.
#'
#' @seealso \code{\link{predict.plsr}}
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#' @noRd
get_validation_metrics <- function(predicted, observed) {
  # predicted  : matrix/data.frame of predictedictions (rows = samples, cols = comps)
  # observed : vector or single-column matrix with ground observed
  P <- as.matrix(predicted)
  Y <- as.matrix(observed)
  
  if (ncol(Y) != 1L)
    stop("`observed` must be a vector or a single-column matrix.")
  
  ok <- complete.cases(Y)
  if (!any(ok))
    return(data.frame(npls = integer(), rmse = numeric(), r2 = numeric()))
  
  P_ok <- P[ok, , drop = FALSE]
  y_ok <- Y[ok, 1, drop = TRUE]
  
  resid <- sweep(P_ok, 1, y_ok, FUN = "-")
  me    <- colMeans(resid)
  rmse  <- sqrt(colMeans(resid^2))
  r2   <- as.vector(cor(P_ok, y_ok)^2)
  
  data.frame(
    ncomp = seq_len(ncol(P)),
    r2   = r2,
    rmse = rmse,
    me   = me,
    row.names = NULL
  )
}

#' @title Tidy and label PLS model components (internal)
#' 
#' @description
#' Internal helper that renames and transposes selected fields of a PLS
#' object to a consistent shape, assigning row and column names based on
#' predictor names and the number of components.
#'
#' @param pls_obj A list-like PLS object containing elements such as
#'   `coefficients`, `projection_mat`, `vip`, `selectivity_ratio`,
#'   `X_loadings`, `Y_loadings`, `weights`, `scores`, and `ncomp`.
#' @param x_names Character vector of predictor variable names.
#'
#' @return The modified `pls_obj` with consistent names and matrix orientations.
#'
#' @details
#' Matrices `coefficients`, `projection_mat`, `vip`, and `selectivity_ratio`
#' are transposed to match downstream conventions (rows = components,
#' cols = variables). Component indices are assigned as `1:ncomp`; scores
#' are labelled with sample and component indices.
#'
#' @seealso [get_plsr()]
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
#' @noRd
assign_pls_names <- function(pls_obj, x_names) {
  # Transpose to convention: rows = components, cols = variables
  pls_obj$coefficients <- t(pls_obj$coefficients)
  pls_obj$projection_mat <- t(pls_obj$projection_mat)
  pls_obj$vip <- t(pls_obj$vip)
  pls_obj$selectivity_ratio <- t(pls_obj$selectivity_ratio)
  
  # Assign column names from predictor names
  colnames(pls_obj$coefficients) <-
    colnames(pls_obj$X_loadings) <-
    colnames(pls_obj$projection_mat) <-
    colnames(pls_obj$vip) <-
    colnames(pls_obj$selectivity_ratio) <-
    colnames(pls_obj$weights) <- x_names
  
  # Assign row names (component indices)
  comp_names <- seq_len(pls_obj$ncomp)
  
  rownames(pls_obj$coefficients) <-
    colnames(pls_obj$bo) <-
    rownames(pls_obj$X_loadings) <-
    rownames(pls_obj$Y_loadings) <-
    rownames(pls_obj$projection_mat) <-
    rownames(pls_obj$vip) <-
    rownames(pls_obj$selectivity_ratio) <-
    rownames(pls_obj$weights) <- comp_names
  
  rownames(pls_obj$bo) <- "Intercepts"
  colnames(pls_obj$Y_loadings) <- "Loadings"
  
  # Scores: rows = samples, cols = components
  colnames(pls_obj$scores) <- comp_names
  rownames(pls_obj$scores) <- seq_len(nrow(pls_obj$scores))
  
  pls_obj
}
