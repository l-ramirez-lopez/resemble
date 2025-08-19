## tools:
## right triangular prism (errror, plsmin, plsmax, k)
##  weighted sd enemble
## vips analysis
## bo analysis
## animation of the vip (reverse from global to local)
## SD vips vs R2 at different neighborhood sizes 

## TODO:
# - add an option to indicate wheter selectivityR and vip must be calculated during building... 
#   if TRUE, opls2 is used (slow) if FALSE opls3 is used (fast)
# 
## downscaling the library to rule-based models


# rgrid <- matrix(NA, max(pls.c), max(pls.c))
# 
# for(pp in 1:105){
#   rgrid[predperformance$minpls[pp],
#         predperformance$maxpls[pp]] <- predperformance$rmse[pp]
# }
# 
# 



# Xr = Xr, 
# Yr = Yr, 
# minF = minF, 
# maxF = maxF, 
# emgrid = emgrid,
# scale = scale, 
# maxiter = maxiter, 
# tol = tol, 
# regression = regression, 
# pc_selection = pc_selection,
# k_kidxmat = kidxmat[ik,],
# k_kidxgrop = kidxgrop[ik,])
# 
# 


#' @title Iterator for Grouped Subsets by Index
#' @description
#' Internal helper that returns an iterator object over subsets of `x` and `y`,
#' grouped by precomputed index sets (e.g., k-nearest neighbors), with optional
#' dissimilarity matrix augmentation.
#'
#' @param x A matrix or data frame of predictors (observations in rows).
#' @param y A vector or matrix of responses, same number of rows as `x`.
#' @param kindx A list or matrix of indices (e.g., kNN indices) used to form 
#' subsets.
#' @param kgroup A vector of group indices into `kindx` to select for each 
#' iteration.
#' @param D (Optional) A square dissimilarity matrix corresponding to `x`, used 
#' to augment the feature set if provided.
#'
#' @return
#' An object of class `isubset`, `abstractiter`, `iter` that yields, on each
#' call to `nextElem()`:
#' \itemize{
#'   \item \code{x}: Subset of predictors (optionally prepended with local 
#'   dissimilarities).
#'   \item \code{y}: Corresponding subset of responses.
#'   \item \code{xval}: The test observation (optionally prepended with 
#'   dissimilarity vector).
#' }
#'
#' @details
#' On each iteration:
#' \itemize{
#'   \item A group of indices is selected using \code{kgroup} from \code{kindx}.
#'   \item Corresponding rows are extracted from \code{x} and \code{y}.
#'   \item If \code{D} is provided, the subset of the dissimilarity matrix for
#'   the current group is used as augmented features.
#' }
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
ith_subsets_by_group <- function(
    x, 
    y, 
    kindx,
    kgroup, 
    D = NULL
){
  
  it_kindx <- iter(kindx, by = "column")
  it_kgroup <- iter(kgroup, by = "column")
  it_xval <- iter(x, by = "row")
  it_D <- iter(D, by = "row")
  
  nextEl <- function() {
    knns <- nextElem(it_kindx)[nextElem(it_kgroup)]
    if(is.null(D)){
      list(x = x[knns, , drop = FALSE], 
           y = y[knns, , drop = FALSE], 
           xval = nextElem(it_xval)) 
    }else{
      idsm <- D[knns, knns]
      list(x = cbind(idsm, x[knns, , drop = FALSE]), 
           y = y[knns, , drop = FALSE], 
           xval = t(c(nextElem(it_D)[knns],
                      nextElem(it_xval))))
    }
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}


#' @title Grouped KNN Subset Iterator with Chunk Size
#'
#' @description
#' Creates an iterator over grouped subsets of rows from `x` and `y` based on
#' nearest neighbor indices (`kindx`) and group selectors (`kgroup`). The
#' iteration is performed in consecutive column batches of size up to
#' `chunk_size` from the input matrices, and optionally incorporates a
#' dissimilarity matrix `D`.
#'
#' @param x A numeric matrix (n × p) representing the full reference dataset.
#' @param y A numeric matrix or vector (n × 1) of associated response values.
#' @param kindx An integer matrix (k × m) where each column contains indices of
#'   the `k` nearest neighbors for the corresponding observation.
#' @param kgroup A logical or integer matrix (k × m) indicating a group
#'   (subset) of rows from `kindx` to be selected for each observation.
#' @param D Optional numeric matrix (n × n) representing a dissimilarity or
#'   distance matrix. If provided, its values for selected neighbors will be
#'   embedded into the `xval` vector.
#' @param chunk_size Integer. The maximum number of columns from `kindx` and
#'   `kgroup` to include in each batch. Defaults to 1. The last batch may be
#'   smaller if the total number of columns is not a multiple of `chunk_size`.
#'
#' @return An iterator object that returns, at each call of `nextElem()`, a list
#'   of lists. Each inner list contains:
#'   - `x`: the reference subset (`k × p` or `k × (p + k)` if `D` used),
#'   - `y`: the corresponding responses (`k × 1`),
#'   - `xval`: the target sample (`1 × p`) or extended with dissimilarities.
#'
#' @details
#' This iterator is useful when local models need to be built for each
#' observation using a selected group of its neighbors (`kgroup`) from `kindx`,
#' and optionally enhanced with pairwise distances via `D`.
#'
#' When `D` is provided, the function will include the distances between the
#' target observation and its selected neighbors in the `xval` element of each
#' result.
#'
#' @seealso [iterators::iter()], [nextElem()]
#'
#' @keywords internal
ith_subsets_by_group_list <- function(
    x, y, kindx, kgroup, D = NULL, chunk_size = 1
) {
  stopifnot(nrow(x) == nrow(y), ncol(kindx) == ncol(kgroup))
  if (!is.null(D)) stopifnot(nrow(D) == nrow(x), ncol(D) == nrow(x))
  
  m <- ncol(kindx)
  chunk_size <- as.integer(chunk_size); stopifnot(chunk_size > 0L)
  
  col_indices <- split(seq_len(m), ceiling(seq_len(m) / chunk_size))
  it_index <- iterators::iter(col_indices)
  
  nextEl <- function() {
    idx <- iterators::nextElem(it_index)
    
    kindx_sub  <- kindx[,  idx, drop = FALSE]
    kgroup_sub <- kgroup[, idx, drop = FALSE]
    
    # Anchor sample per column = first row of kindx
    anchor_idx <- as.integer(kindx_sub[1L, ])
    xval_sub   <- x[anchor_idx, , drop = FALSE]
    
    chunks <- vector("list", length(idx))
    for (i in seq_along(idx)) {
      sel <- as.logical(kgroup_sub[, i]); sel[is.na(sel)] <- FALSE
      sel[1L] <- FALSE  # never include anchor among neighbours
      
      knns <- as.integer(kindx_sub[sel, i])
      
      if (is.null(D)) {
        chunks[[i]] <- list(
          x = x[knns, , drop = FALSE],
          y = y[knns, , drop = FALSE],
          xval = xval_sub[i, , drop = FALSE]
        )
      } else {
        idsm  <- D[knns, knns, drop = FALSE]
        xdata <- cbind(idsm, x[knns, , drop = FALSE])
        chunks[[i]] <- list(
          x = xdata,
          y = y[knns, , drop = FALSE],
          xval = t(c(D[anchor_idx[i], knns], xval_sub[i, , drop = FALSE]))
        )
      }
    }
    chunks
  }
  
  obj <- list(nextElem = nextEl)
  class(obj) <- c("isubsetgroup", "abstractiter", "iter")
  obj
}




#' @title Iterator for Subsets by Index
#' @description
#' Internal helper to iterate over subsets of `x` and `y` based on `kindx`
#' indices. Optionally augments `x` with a dissimilarity matrix `D` and
#' filters out rows where `y` is NA.
#'
#' @param x A matrix or data frame of predictors (rows = observations).
#' @param y A vector or matrix of responses (same rows as `x`).
#' @param kindx A list or matrix of row indices (e.g., nearest neighbors).
#' @param D Optional square dissimilarity matrix matching rows of `x`.
#'
#' @return
#' An iterator object of class `isubset`, `abstractiter`, `iter`. Each
#' call to `nextElem()` yields:
#' \itemize{
#'   \item \code{x}: A subset of `x`, optionally prepended with dissimilarities.
#'   \item \code{y}: Corresponding subset of `y`, with NA rows removed.
#' }
#'
#' @details
#' For each iteration, indices from `kindx` are used to extract rows
#' from `x` and `y`. If `D` is given, `x` is augmented with local
#' dissimilarities from `D[knns, knns]`. Observations where `y` is NA
#' are excluded from both `x` and `y` in the returned subset.
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
ith_subsets <- function(
    x, 
    y, 
    kindx,
    D = NULL
){
  it_kindx <- iter(kindx, by = "column")
  it_xval <- iter(x, by = "row")
  it_D <- iter(D, by = "row")
  
  nextEl <- function() {
    knns <- nextElem(it_kindx)
    ixr <- x[knns, ]
    iyr <- y[knns, , drop = FALSE]
    if (!is.null(D)) {
      ixr <- cbind(D[knns, knns], ixr)
    }
    list(
      x = ixr[!is.na(iyr), ], 
      y = iyr[!is.na(iyr), , drop = FALSE]
    )
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}



#' @title Chunked Iterator for Subsetting by Nearest Neighbors
#' @description
#' Creates an iterator that returns batched subsets of observations for
#' local modeling. Each batch corresponds to a set of consecutive indices
#' (rows) of size up to `chunk_size`. For each observation in the batch,
#' its corresponding nearest neighbors are used to extract the `x`, `y`,
#' and optional `D` values.
#'
#' @param x A numeric matrix of predictors (n × p), where `n` is the number of
#'   observations and `p` the number of features.
#' @param y A numeric vector or matrix of responses of length `n` or
#'   dimension `n × 1`.
#' @param kindx An integer matrix of size `k × n`, where each column contains
#'   the indices of the `k` nearest neighbors for the corresponding observation.
#' @param D Optional numeric matrix of pairwise dissimilarities (n × n).
#'   If provided, for each subset the dissimilarities among the selected
#'   neighbors are prepended as additional columns to the predictors.
#' @param chunk_size Integer. Maximum number of observations per batch.
#'   The `n` observations are split into consecutive batches of size
#'   `chunk_size`, except for the final batch which may be smaller.
#'
#' @return An iterator object compatible with `nextElem()` from the
#'   **iterators** package. Each call to `nextElem()` returns a list of
#'   sublists (one per observation in the batch), where each sublist contains:
#'   - `x`: predictor matrix for that local subset (with optional dissimilarity)
#'   - `y`: response vector/matrix for that local subset
#'
#' @examples
#' \dontrun{
#' library(iterators)
#' it <- ith_subsets_list(Xr, Yr, kindx, D, chunk_size = 100)
#' 
#' while (TRUE) {
#'   chunk <- tryCatch(nextElem(it), error = function(e) {
#'     if (inherits(e, "StopIteration")) return(NULL)
#'     stop(e)
#'   })
#'   if (is.null(chunk)) break
#'   print(length(chunk))  # Number of subsets in the batch
#' }
#' }
#'
#' @keywords internal
ith_subsets_list <- function(x, y, kindx, D = NULL, chunk_size = 1) {
  stopifnot(nrow(x) == nrow(y))
  if (!is.null(D)) stopifnot(nrow(D) == nrow(x), ncol(D) == nrow(x))
  
  # We iterate over columns of kindx (one local model per column)
  m <- ncol(kindx)
  stopifnot(is.numeric(chunk_size), length(chunk_size) == 1L)
  chunk_size <- as.integer(chunk_size)
  stopifnot(chunk_size > 0L)
  
  # Build consecutive column-index chunks of size `chunk_size`
  col_chunks <- split(seq_len(m), ceiling(seq_len(m) / chunk_size))
  it_chunks <- iterators::iter(col_chunks)
  
  nextEl <- function() {
    idx <- iterators::nextElem(it_chunks)  # indices of columns in `kindx`
    
    lapply(idx, function(i) {
      knns <- as.integer(kindx[, i])
      ixr  <- x[knns, , drop = FALSE]
      iyr  <- y[knns, , drop = FALSE]
      
      if (!is.null(D)) {
        ixr <- cbind(D[knns, knns, drop = FALSE], ixr)
      }
      
      keep <- !is.na(iyr)
      list(
        x = ixr[keep, , drop = FALSE],
        y = iyr[keep, , drop = FALSE]
      )
    })
  }
  
  obj <- list(nextElem = nextEl)
  class(obj) <- c("chunked_isubset", "abstractiter", "iter")
  obj
}




#' @title Iterator for Prediction Subsets (PLS Local Models)
#' @description
#' Internal helper to iterate over prediction inputs for local PLS models.
#'
#' @param plslib Matrix of PLS regression vectors (one per reference spectrum).
#' @param Xu Matrix of new observations to predict (rows = samples).
#' @param xunn Matrix of nearest neighbor indices for each row of Xu.
#' @param xscale Matrix of scaling vectors for each reference model.
#' @param dxrxu Optional dissimilarity matrix (ref vs. new samples).
#'
#' @return
#' An iterator of class `isubset`, `abstractiter`, `iter`. Each call
#' to `nextElem()` yields:
#' \itemize{
#'   \item \code{iplslib}: Local PLS regression vectors for neighbors.
#'   \item \code{ixscale}: Scaling vectors for those models.
#'   \item \code{ixunn}: Indices of nearest neighbors for prediction.
#'   \item \code{ixu}: New observation from `Xu` to predict.
#'   \item \code{idxrxu}: (Optional) dissimilarity vector to neighbors.
#' }
#'
#' @details
#' For each new sample in `Xu`, retrieves its nearest neighbors, along
#' with corresponding PLS vectors, scaling vectors, and optional
#' dissimilarity values. Intended for local PLS prediction routines.
#' @author Leonardo Ramirez-Lopez
#' @keywords internal
ith_pred_subsets <- function(
    plslib,
    Xu,
    xunn,
    xscale,
    dxrxu = NULL
){
  it_xunn <- iter(xunn, by = "column")
  it_xu <- iter(Xu, by = "row")
  
  if (!is.null(dxrxu)) {
    it_dxrxu <- iter(dxrxu, by = "column")
  } else {
    it_dxrxu <- NULL
  }
  
  nextEl <- function() {
    ixunn <- nextElem(it_xunn)
    ixu <- nextElem(it_xu)
    
    if (!is.null(it_dxrxu)) {
      idxrxu <- nextElem(it_dxrxu)
      idxrxu <- idxrxu[ixunn]
    } else {
      idxrxu <- NULL
    }

    ixscale <- xscale[ixunn, ]
    iplslib <- plslib[ixunn, ]
    
    list(
      iplslib = iplslib,
      ixscale = ixscale,
      ixunn = ixunn,
      ixu = ixu,
      idxrxu = idxrxu
    )
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}

# knns <- kindx[,..ii..]
# iyr <- y[knns, ,drop = FALSE]
# ixr <- cbind(D[knns, knns], x[knns, ])
# ixr <- ixr[!is.na(iyr), ]
# iyr <- iyr[!is.na(iyr), , drop = FALSE]
# 




#' #' @title Quantile Stats of Neighbor Responses
#' #' @description
#' #' Computes quantiles of the response values among the filtered
#' #' nearest neighbors for a specific sample.
#' #'
#' #' @param ..i.. Index of the sample being analyzed.
#' #' @param kidxmat Matrix of nearest neighbor indices (rows = neighbors).
#' #' @param kidxgrop Matrix mask to filter out same-group neighbors (TRUE = keep).
#' #' @param Yr Numeric response vector or 1-col matrix.
#' #' @param k Vector of neighbor counts for each sample.
#' #'
#' #' @return
#' #' A matrix of quantiles (0%, 5%, 25%, 50%, 75%, 95%, 100%) of the
#' #' response values among the valid nearest neighbors.
#' #'
#' #' @details
#' #' For the i-th sample, this function selects the top-k neighbors for
#' #' each sample, filters using the group mask, and computes quantile
#' #' statistics from their corresponding response values in `Yr`.
#' #' @author Leonardo Ramirez-Lopez
#' #' @keywords internal
#' i_nn_stats <- function(..i.., kidxmat, kidxgrop, Yr, k) {
#'   ik <- k[..i..]
#'   
#'   jstats <- sapply(
#'     1:ncol(kidxmat),
#'     FUN = function(..j.., kidxmat, kidxgrop, Yr, ik){
#'       inn <- kidxmat[1:ik,..j..]
#'       inn <- inn[kidxgrop[1:ik,..j..]]
#'       quantile(Yr[inn,], c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.00))
#'       #sd(Yr[inn,])
#'     }, 
#'     kidxgrop = kidxgrop,
#'     kidxmat = kidxmat,
#'     Yr = Yr,
#'     ik = ik
#'   )
#'   t(jstats)
#' }

#' @title Internal: Fit weighted local PLS models for all reference samples
#'
#' @description
#' Computes predictions from weighted local PLS models for each observation in
#' the reference set `Xr`, based on its nearest neighbors and optionally using
#' a dissimilarity matrix as extra predictive information.
#'
#' @param ..k.. Index (integer) indicating the number of neighbors to use
#'              from vector `k`.
#' @param Xr A matrix or data.frame of predictor variables for the reference
#'           data (rows = observations, columns = variables).
#' @param Yr A numeric vector of response values corresponding to `Xr`.
#' @param k A vector indicating the number of nearest neighbors to test for
#'          each reference observation.
#' @param min_component Minimum number of components to use in PLS regression.
#' @param max_component Maximum number of components to use in PLS regression.
#' @param emgrid A matrix or grid of external evaluation weights.
#' @param scale Logical; if TRUE, scale predictors before regression.
#' @param maxiter Maximum number of iterations for the PLS algorithm.
#' @param tol Convergence tolerance for the PLS algorithm.
#' @param pc_selection A list specifying the principal component selection
#'        method. Must contain:
#' \itemize{
#'   \item{\code{method}:}{Character string; one of \code{"opc"},
#'         \code{"cumvar"}, \code{"var"}, or \code{"manual"}.}
#'   \item{\code{value}:}{Numeric; threshold or fixed number of components,
#'         depending on the method.}
#' }
#' For example, \code{list(method = "opc", value = 40)}.
#'
#' @param kidxmat Matrix of neighbor indices for each observation (columns =
#'        samples, rows = neighbors).
#' @param kidxgrop Logical matrix indicating which neighbors belong to the
#'        same group as the target and should be excluded.
#' @param dissimilarity_mat Optional matrix of dissimilarities between
#'        reference observations, used to enhance feature input for PLS.
#' @param pb A progress bar object (e.g., from \code{txtProgressBar()}).
#' @param ... Further arguments (currently unused).
#'
#' @return A matrix of predicted values (columns = observations in `Xr`,
#'         rows = emgrid dimensions).
#'
#' @details
#' For each reference observation in `Xr`, the function builds a local
#' training subset using the `ith_subsets_by_group()` function, fits a local
#' weighted PLS model via `ith_local_fit()`, and produces predictions based on
#' the evaluation grid (`emgrid`). Dissimilarity values can be optionally used
#' to augment predictive inputs. Parallel computation is handled via
#' `foreach` + `%dopar%`.
#'
#' @author Leonardo Ramirez-Lopez
#'
#' @keywords internal
.get_all_fits <- function(
    ..k.., 
    Xr, 
    Yr, 
    k,
    min_component, 
    max_component, 
    emgrid,
    scale, 
    maxiter, 
    tol, 
    pc_selection,
    kidxmat,
    kidxgrop, 
    dissimilarity_mat = NULL,
    pb,
    chunk_size,
    ...
) {
  ik <- seq_len(k[..k..])
  
  setTxtProgressBar(pb, ..k..)
  
  # ---- Optionally setup parallel backend ----
  # n.cores <- 20
  # clust <- parallel::makeCluster(n.cores)
  # doParallel::registerDoParallel(clust)
  
  # ith_preds_template <- matrix(
  #   NA, 
  #   nrow = nrow(emgrid), 
  #   ncol = ceiling(nrow(Xr) / n_chunks)
  # )
  # ---- Parallel foreach execution (default mode) ----
  
  n_chunks <- nrow(Xr) / chunk_size 
  
  n_xr <- nrow(Xr)
  ith_preds_template <- matrix(
    NA,
    nrow = nrow(emgrid),
    ncol = chunk_size 
  )
  innpreds <- foreach(
    i = 1:nrow(Xr), 
    ksubsets = ith_subsets_by_group_list(
      x = Xr, 
      y = Yr, 
      kindx = kidxmat[ik, ], 
      kgroup = kidxgrop[ik, ], 
      D = dissimilarity_mat, 
      chunk_size = chunk_size
    ),
    .export = c(
      "ith_pred_subsets",
      "ith_subsets",
      "ith_subsets_by_group",
      ".get_all_fits",
      "ith_local_fit",
      "opls"
    )
  ) %dopar% {
    ith_preds <- ith_preds_template
    for (j in 1:length(ksubsets)) {
      ith_preds[, j] <- ith_local_fit(
        X = ksubsets[[j]]$x,
        Y = ksubsets[[j]]$y,
        xval = ksubsets[[j]]$xval,
        emgrid = emgrid,
        max_component = max_component,
        min_component = min_component,
        scale = scale,
        maxiter = maxiter,
        tol = tol
      )
    }
    if (j < ncol(ith_preds)) {
      ith_preds <- ith_preds[, 1:j, drop = FALSE]
    }
    ith_preds
  }
  innpreds <- do.call("cbind", innpreds)
  return(innpreds)
}

# #' @title Internal: Fit a local weighted PLS model and predict for a query point
# #'
# #' @description
# #' Fits a local Partial Least Squares (PLS) model using a neighborhood subset
# #' and computes a weighted prediction for a target sample. The weighting is
# #' done over multiple components using a provided evaluation grid.
# #'
# #' @param ilocalsubset A list with elements:
# #' \itemize{
# #'   \item{\code{x}}: matrix of predictors from the local neighborhood.
# #'   \item{\code{y}}: vector or 1-col matrix of corresponding responses.
# #'   \item{\code{xval}}: query sample (usually one row) to predict.
# #' }
# #' @param min_component Minimum number of PLS components to use in prediction.
# #' @param max_component Maximum number of PLS components to fit.
# #' @param emgrid A numeric matrix used to weight component-wise predictions.
# #' @param scale Logical; whether to scale predictors before PLS fitting.
# #' @param maxiter Maximum number of iterations for the PLS algorithm.
# #' @param tol Convergence tolerance for the PLS algorithm.
# #' @param pc_selection A list defining the principal component selection
# #'        strategy. (Included for consistency, not used directly here.)
# #' @param ... Additional arguments (currently unused).
# #'
# #' @return
# #' A numeric matrix of weighted predictions (rows = grid rows, 1 column).
# #'
# #' @details
# #' The function performs the following steps:
# #' \enumerate{
# #'   \item Fits a PLS model on `ilocalsubset$x` and `ilocalsubset$y` with
# #'         up to `max_component` components.
# #'   \item Computes component weights for `ilocalsubset$xval` using
# #'         \code{get_local_pls_weights()}.
# #'   \item Predicts component-wise responses using \code{predict_opls()}.
# #'   \item Applies a weighted average of predictions using `emgrid` and the
# #'         component weights.
# #' }
# #' The output is a matrix of weighted predictions using the weighted model grid.
# #'
# #' @author Leonardo Ramirez-Lopez
# #'
# #' @keywords internal
# ith_local_fit <- function(
    #     ilocalsubset,
#     min_component,
#     max_component,
#     emgrid,
#     scale, 
#     maxiter, 
#     tol, 
#     pc_selection,
#     ...
# ){ 
#   
#   ipreds <- run_ipls_prediction(
#     X = ilocalsubset$x,
#     Y = ilocalsubset$y,
#     xval = ilocalsubset$xval,
#     emgrid = emgrid,
#     max_component = max_component,
#     min_component = min_component,
#     scale = scale,
#     maxiter = maxiter,
#     tol = tol
#   )
#   ipreds
# }


# #' @title Internal: Compute weighted final PLS model outputs
# #' 
# #' @description
# #' Fits a PLS model on a local neighborhood and computes a weighted combination
# #' of regression coefficients, intercept, VIP scores, selectivity ratios, and
# #' predictor scales based on component-specific weights.
# #' 
# #' @param ilocalsubset A list with local subset data:
# #' \itemize{
# #'   \item \code{x}: matrix of predictors from the local neighborhood.
# #'   \item \code{y}: response vector or 1-col matrix for local samples.
# #' }
# #' @param min_component Integer; minimum number of components used for output.
# #' @param max_component Integer; maximum number of components to fit in PLS.
# #' @param scale Logical; whether to scale predictors before PLS fitting.
# #' @param maxiter Integer; maximum number of iterations for the PLS algorithm.
# #' @param tol Numeric; convergence tolerance for the PLS algorithm.
# #' @param ... Additional arguments (currently unused).
# #' 
# #' @return
# #' A numeric vector containing the following, in order:
# #' \itemize{
# #'   \item Intercept term (weighted average of component intercepts)
# #'   \item Weighted regression coefficients
# #'   \item Weighted VIP scores (normalized by column SDs)
# #'   \item Weighted selectivity ratios (normalized by column SDs)
# #'   \item Predictor scaling values (from training transformation)
# #' }
# #' 
# #' @details
# #' This function fits a PLS model on the local neighborhood using
# #' \code{opls_get_all()}, computes local weights with
# #' \code{get_local_pls_weights()}, and combines all component-wise metrics
# #' into a final weighted output.
# #' 
# #' Component indices are selected using the range
# #' \code{min_component:max_component}. Weighting is applied to each component
# #' across coefficients, VIPs, and selectivity ratios. All component-wise scores
# #' are normalized before aggregation. The result can be used as a compressed
# #' representation of the final prediction model.
# #' 
# #' Note: the prediction target is assumed to be \code{x[1, ]} of the subset.
# #' 
# #' @author Leonardo Ramirez-Lopez
# #' 
# #' @keywords internal
# final_fits <- function(
    #     ilocalsubset,
#     min_component, 
#     max_component, 
#     scale, 
#     maxiter, 
#     tol, 
#     ...
# ){
#   
#   result <- final_fits_cpp(
#     X = ilocalsubset$x,
#     Y = ilocalsubset$y,
#     new_x = ilocalsubset$x[1, , drop = FALSE],
#     min_component,
#     max_component,
#     scale = scale,
#     maxiter = maxiter,
#     tol = tol
#   )
#   unlist(result, recursive = FALSE, use.names = FALSE)
# }


#' @title Internal: Predict using local PLS coefficients and optional
#' dissimilarities
#'
#' @description
#' Computes a prediction for a new observation using a local PLS model
#' represented by coefficients (`plslib`). The prediction is based on
#' inverse-scaled feature values. If a dissimilarity vector is provided, it is
#' appended to the input features before inverse scaling.
#'
#' @param plslib A matrix of PLS model coefficients. First column is the
#'        intercept; remaining columns are coefficients for scaled features.
#' @param xscale A matrix of scaling values for each feature (same dimensions
#'        as \code{plslib[,-1]}).
#' @param Xu A numeric vector (or single-row matrix) representing the query
#'        sample to be predicted.
#' @param xunn Currently unused; placeholder for nearest neighbor indices.
#' @param dxrxu Optional numeric vector of dissimilarities between `Xu` and
#'        reference samples. If provided, used as additional predictive
#'        features.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric scalar containing the predicted response for `Xu`.
#'
#' @details
#' The function constructs a synthetic predictor vector by dividing
#' `xscale` by `Xu` (or by the concatenation of `dxrxu` and `Xu`), then
#' inverting the result. This reflects a distance- or similarity-based
#' transformation. The resulting vector is then multiplied element-wise with
#' the PLS coefficients, and the sum is added to the intercept to compute
#' the final prediction.
#'
#' @author Leonardo Ramirez-Lopez
#'
#' @keywords internal
ith_pred <- function(plslib, xscale, Xu, xunn, dxrxu = NULL, ...){
  if (is.null(dxrxu)) {
    ixus <- sweep(xscale, MARGIN = 2, STATS =  Xu, FUN = "/" )
    ixus <- 1 / ixus
    ipred <- plslib[, 1] + rowSums(plslib[, -1] * ixus)
  } else {
    dxu <- c(dxrxu, Xu)
    ixus <- sweep(xscale, MARGIN = 2, STATS =  dxu, FUN = "/" )
    ixus <- 1 / ixus
    ipred <- plslib[, 1] + rowSums(plslib[, -1] * ixus)
  }
  ipred
}




#########################################################################

# 
# diss.predictors = FALSE
# 
# 
# 
# ## new args
# k <- seq(40, 150, 20)
# group <- as.factor(paste("g", 1:nrow(Xr), sep = ""));# group[c(1,2,4)] <- "g1"
# pls.c <- c(4, 17)
# pls_max_iter = 1
# pls_tol = 1e-6
# metric = "rmse"
# 
# data.frame(ix = 1:length(group), group)




#' Predict from a funlib object
#'
#' Generate predictions using a local ensemble model from a \code{funlib} object.  
#' For each new sample, a neighbourhood of local models is retrieved based on  
#' dissimilarity, and predictions are obtained via a weighted combination of  
#' these models.
#'
#' @param object A fitted object of class \code{"funlib"}, typically created by
#'   a function that builds a function library of local models.
#' @param newdata A data frame or matrix containing new predictor values.
#'   Must include all predictors required by \code{object}.
#' @param weighting Character string specifying the kernel weighting function 
#'   applied to neighbours when combining predictions. Options are: 
#'   \code{"triweight"} (default), \code{"tricube"}, \code{"triangular"}, 
#'   \code{"quartic"}, \code{"parabolic"}, or \code{"gaussian"}.
#' @param probs A numeric vector of probabilities in \[0, 1\]. 
#'   (Currently not implemented.) Intended to specify weighted quantiles of the  
#'   predictions across models to provide an empirical approximation to  
#'   prediction intervals. Default is \code{seq(0, 1, 0.25)}.
#' @param range_pred_lim Logical. If \code{TRUE}, predictions falling outside
#'   the 5–95% range of neighbour target values are clipped to those limits.
#' @param local Logical; if \code{TRUE}, predictions are based on local models.
#' @param residual_cutoff Numeric. Threshold for excluding models with residuals
#'   greater than this value. Default is \code{Inf} (no exclusion).
#' @param diss_method Optional character string specifying the dissimilarity 
#'   method to use when retrieving models. Valid options include 
#'   \code{"pca"}, \code{"pca.nipals"}, \code{"pls"}, \code{"mpls"}, 
#'   \code{"cor"}, \code{"euclid"}, \code{"cosine"}, or \code{"sid"}.
#' @param ... Further arguments passed to the dissimilarity computation.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{predictions}}{A data frame with the following columns:
#'     \describe{
#'       \item{\code{pred}}{Predicted values}
#'       \item{\code{pred_sdev}}{Weighted standard deviation of predictions}
#'       \item{\code{q[probability]}}{Weighted quantiles of predictions}
#'       \item{\code{gh}}{(if available) Global heterogeneity dissimilarity score}
#'       \item{\code{min_yr_in_neighborhood}}{Minimum neighbour response (5th percentile)}
#'       \item{\code{max_yr_in_neighborhood}}{Maximum neighbour response (95th percentile)}
#'       \item{\code{below_lower_lim}, \code{above_upper_lim}}{Logicals indicating 
#'         whether the prediction falls outside these limits}
#'     }
#'   }
#'   \item{\code{neighbour_list}}{A list with components:
#'     \describe{
#'       \item{\code{indices}}{Indices of the selected neighbours}
#'       \item{\code{diss_scores}}{Their dissimilarity scores}
#'     }
#'   }
#' }
#' @details
#' This function applies dissimilarity-based local modelling to compute 
#' ensemble predictions for new data. For each query sample, a neighbourhood of 
#' \code{optimalk} local models is selected using the specified dissimilarity 
#' method. Predictions from these models are then combined using the selected 
#' kernel weighting function.
#'
#' @references
#' Cleveland, W. S., & Devlin, S. J. (1988). Locally weighted regression: 
#' An approach to regression analysis by local fitting. 
#' \emph{Journal of the American Statistical Association}, 83(403), 596–610.
#'
#' Næs, T., & Isaksson, T. (1990). Locally weighted regression and scatter 
#' correction methods in practical spectroscopy. \emph{Applied Spectroscopy}, 
#' 44(3), 378–383.
#'
#' @author Leonardo Ramirez-Lopez
#' @export

predict.funlib <- function(
    object, 
    newdata, 
    weighting = "triweight",
    probs = seq(0, 1, 0.25),
    range_pred_lim = FALSE, 
    local = TRUE, 
    residual_cutoff = Inf,
    diss_method = NULL, 
    ...
) {
  
  if (object$dissimilatity$diss_method == "Precomputed dissimilarity matrix" & is.null(diss_method)) {
    stop(
      "The models were built using a precomputed dissimilarity matrix. ",
      "To select models, you must explicitly specify a dissimilarity method ",
      "in the 'diss_method' argument. Valid options include: ",
      "'pca', 'pca.nipals', 'pls', 'mpls', ",
      "'cor', 'euclid', 'cosine', or 'sid'."
    )
  } else {
    diss_method_type <- object$dissimilatity$diss_method
  }
  
  if (!is.null(diss_method)) {
    diss_method_type <- match.arg(
      diss_method, 
      c("pca", "pca.nipals", "pls", "mpls", "cor", "euclid", "cosine", "sid")
    )
  }

  if(!"funlib" %in% class(object)){
    stop("object must be of class 'funlib'")
  }
  
  if(!all(colnames(object$functionlibrary$B) %in% colnames(newdata))){
    stop("missing predictor variables in newdata")
  }
  
  newdata <- newdata[, colnames(newdata) %in% colnames(object$functionlibrary$B)]

  ghd <- NULL
  if (diss_method_type %in% c("pca", "pls")) {
    scnew <- predict(object$dissimilatity$projection, newdata)
    zcenter <- resemble:::get_column_means(object$dissimilatity$projection$scores) 
    zscale <- resemble:::get_column_sds(object$dissimilatity$projection$scores)
    
    dsmxu <- dissimilarity(
      Xr = scale(
        object$dissimilatity$projection$scores[object$anchor_indices, ], 
        center = zcenter, 
        scale = zscale
      ),
      Xu = scale(
        scnew, 
        center = zcenter, 
        scale = zscale
      ),
      diss_method = "euclid",
      center = FALSE,
      scaled = FALSE,
      gh = FALSE
    )
    
  } else {

    dsmxu <- dissimilarity(
      Xr = scale(
        object$scale$local.x.center, 
        center = object$scale$centre, scale = object$scale$scale
      ),
      Xu = scale(
        newdata, 
        center = object$scale$centre, 
        scale = object$scale$scale
      ),
      diss_method = diss_method_type,
      center = FALSE,
      scaled = FALSE,
      gh = FALSE,
      ...
    )
    
  }
  
  if (!is.null(object$gh)) {
    if (diss_method_type != "pls") {
      scnew <- predict(object$gh$projection, newdata)
    }
    ghd <- dissimilarity(
      Xr = scnew,
      Xu = t(colMeans(object$gh$projection$scores)),
      diss_method = "euclid",
      center = TRUE,
      scaled = TRUE,
      gh = FALSE
    )$dissimilarity
    ghd <- as.vector(ghd)
  }
  
  dots <- list(...)
  if (diss_method_type == "cor") {
    if ("ws" %in% names(dots)) {
      cat(
        "Retriving models using correlation dissimilarity with a window size of ", ws, "...\n"
      )
    } else {
      cat("Retriving models using correlation dissimilarity with full window...\n")
    }
  } else {
    cat("Retriving models using ", diss_method_type, " dissimilarity...\n")
  }
  # Weights are defined according to a tricubic function 
  # as in Cleveland and Devlin (1988) and Naes and Isaksson (1990).
  xunn <- apply(
    dsmxu$dissimilarity, MARGIN = 2, FUN = order
  )[1:object$sel.param$optimalk,]
  
  xudss <- apply(
    dsmxu$dissimilarity, 
    MARGIN = 2, 
    FUN = sort
  )[1:object$sel.param$optimalk, ]
  
  
  ## this might become an argument to cancel models with high residuals (rd > xx)
  rd <- residual_cutoff
  ## this is for cancelling models with high residuals (rd > xx)
  if (!is.null(object$residuals)) {
    res <- object$residuals
    # res <- abs(scale(res, center = TRUE, scale = TRUE))
    plot(res)
    res[res < rd] <- 0
    res[res >= rd] <- 1
    res <- as.logical(res)
  } else {
    res <- rep(FALSE, ncol(xunn))
  }
  ## this deactivates model cancelling (for the moment)
  # res[] <- FALSE 
  
  xudss2 <- sapply(
    1:ncol(xunn),
    FUN = function(diss, nn, cs, ..i..){
      d1 <- diss[nn[,..i..],..i..]
      d2 <- cs[nn[,..i..]]
      d1[d2] <- max(d1)
      d1
    },
    diss = dsmxu$dissimilarity, 
    nn = xunn,
    cs = res
  )[1:object$sel.param$optimalk, ]
  
  dweights <- sweep(
    xudss2, 
    MARGIN = 2, 
    STATS = get_column_maxs(xudss2), 
    FUN = "/", 
    check.margin = FALSE
  )
  
  
  
  if (weighting == "triweight") {
    # Triweight 
    dweights <- (1 - dweights^2)^3
    #curve((1 - x^2)^3, from = 0, to = 1, ylab = "parabolic")
  }
  
  if (weighting == "tricube") {
    dweights <- (1 - dweights^3)^3
    #curve((1 - x^3)^3, from = 0, to = 1, ylab = "parabolic")
  }
  
  if (weighting == "triangular"){
    # Triweight 
    dweights <- (1 - dweights)
    #curve((1 - x), from = 0, to = 1, ylab = "parabolic")
  }
  
  if (weighting == "quartic") {
    dweights <- (1 - dweights^2)^2
    #curve((1 - x^2)^2, from = 0, to = 1, ylab = "parabolic")
  }
  
  if (weighting == "parabolic") { 
    dweights <- (1 - dweights^2)
    #curve((1 - x^2),from = 0, to = 1, ylab = "parabolic")
  }
  
  if (weighting == "gaussian") {
    dweights <- exp(-dweights^2)
    #curve(exp(-x^2), from = 0, to = 1, ylab = "gaussian")
  }
  
  dweights <- sweep(
    dweights, 
    MARGIN = 2, STATS = colSums(dweights), FUN = "/", 
    check.margin = FALSE
  )
  
  
  if ("Bk" %in% names(object$functionlibrary)) {
    plslib <- cbind(
      object$functionlibrary$B0,
      object$functionlibrary$Bk,
      object$functionlibrary$B
    )
    localscale <- cbind(object$scale$local.diss.scale, object$scale$local.x.scale)
    xudiss <- dsmxu$dissimilarity
  } else {
    plslib <- cbind(
      object$functionlibrary$B0,
      object$functionlibrary$B
    )
    localscale <- object$scale$local.x.scale
    xudiss <- NULL
  }
  
  localset <- ith_pred_subsets(
    plslib = plslib,
    Xu = newdata,
    xunn = xunn,
    xscale = localscale,
    dxrxu = xudiss
  )
  
  cat("Computing predictions...\n")
  yupreds <- foreach(
    i = 1:nrow(newdata), 
    .export = c(
      "ith_pred_subsets",
      "ith_subsets",
      "ith_subsets_by_group",
      ".get_all_fits",
      "ith_local_fit",
      "ith_pred"
    ),
    iset = localset
  ) %dopar% { 
    ipd <- ith_pred(
      plslib = iset$iplslib, 
      xscale = iset$ixscale, 
      Xu = iset$ixu, 
      xunn = iset$ixunn,
      dxrxu = iset$idxrxu
    )
    ipd
  }
  yupreds <- do.call("rbind", yupreds)
  
  ## weighted standard deviation
  # xixm <- sweep(
  #   yupreds, 
  #   MARGIN = 1, 
  #   STATS = rowMeans(yupreds), 
  #   FUN = "-", 
  #   check.margin = FALSE
  # )^2
  # a <- rowSums(xixm * t(dweights))
  # b <- (colSums(!dweights == 0) - 1) / colSums(!dweights == 0)
  # weighted_sdev <- (a / b)^0.5

  # weighted mean per row
  dweights <- t(dweights)
  weighted_yu_preds <- yupreds * dweights

  wpredictions <- rowSums(weighted_yu_preds)
  
  # Centered squared deviations (yupreds already weighted)
  sq_dev <- (weighted_yu_preds - wpredictions)^2
  
  # Apply same weights again
  wvar <- rowSums(sq_dev)
  
  # Standard deviation = sqrt(variance)
  weighted_sdev <- sqrt(wvar)
  wquantiles <- weighted_quantiles(yupreds[, -ncol(dweights)], w = dweights[, -ncol(dweights)], probs = probs)
    
  preds <- data.frame(
    pred = rowSums(weighted_yu_preds), 
    pred_sdev = weighted_sdev, 
    pred_quantiles = wquantiles
  )
  
  preds$gh <- ghd
  
  if (!is.null(rownames(newdata))) {
    rownames(preds) <- rownames(newdata)
  }
  
  
  preds <- list(predictions = preds)
  preds$neighbour_list <- list(
    indices = xunn, 
    diss_scores = xudss
  )
  
  stst <- object$yu.nnstats[[paste("k.",object$sel.param$optimalk, sep = "")]]
  
  liminf <- sapply(1:nrow(newdata), FUN = function(x, ni, ..i..) min(x[ni[,..i..], "5%"]), x = stst, ni = xunn)
  limsup <- sapply(1:nrow(newdata), FUN = function(x, ni, ..i..) max(x[ni[,..i..], "95%"]), x = stst, ni = xunn)
  
  preds$predictions$min_yr_in_neighborhood <- liminf
  preds$predictions$max_yr_in_neighborhood <- limsup
  
  preds$predictions$below_lower_lim <- preds$predictions$pred < preds$predictions$min_yr_in_neighborhood
  preds$predictions$above_upper_lim <- preds$predictions$pred > preds$predictions$max_yr_in_neighborhood
  
  if (range_pred_lim) {
    preds$predictions$pred[preds$predictions$below_lower_lim] <- preds$predictions$min_yr_in_neighborhood[preds$predictions$below_lower_lim]
    preds$predictions$pred[preds$predictions$above_upper_lim] <- preds$predictions$max_yr_in_neighborhood[preds$predictions$above_upper_lim]
  }
  
  return(preds)
}


# predict.validationfunlib <- function(object,...){
#   stop("this object only contains validation data")
# }



.get_coefficient_sds <- function(object){
  if(!"funlib" %in% class(object)){
    stop("object is not of class 'funlib'")
  }
  nf <- ifelse(
    is.null(object$functionlibrary$Bk), 
    ncol(object$functionlibrary$B), 
    ncol(object$functionlibrary$B) + ncol(object$functionlibrary$Bk)
  )
  
  b0 <- object$functionlibrary$B0/nf 
  
  stdb <- NULL
  stdb$b <- sweep(object$functionlibrary$B,
                  MARGIN = 1,
                  FUN = "+",
                  STATS = b0)
  
  
  if(!is.null(object$functionlibrary$Bk)){
    stdb$bk <- sweep(object$functionlibrary$Bk,
                     MARGIN = 1,
                     FUN = "+",
                     STATS = b0)
  }  
  return(stdb)
}


#' Compute weighted quantiles
#'
#' Computes weighted quantiles from a numeric vector and associated weights.
#' Observations are sorted, weights are normalized to sum to 1, and
#' cumulative weights are used to determine quantile positions via linear
#' interpolation.
#'
#' @param x A numeric vector of observations.
#' @param w A numeric vector of non-negative weights, the same length as `x`.
#'   Weights are internally normalized to sum to 1.
#' @param probs A numeric vector of probabilities in [0, 1] specifying the
#'   quantile levels to compute (e.g., `c(0.25, 0.5, 0.75)`).
#'
#' @return A numeric vector of weighted quantiles corresponding to `probs`.
#'
#' @details
#' The function uses a weighted version of the empirical cumulative
#' distribution function (CDF). It sorts `x` and applies normalized weights
#' to compute cumulative sums, which are then used to interpolate the
#' quantiles at specified `probs`. Ties in cumulative weights are averaged.
#'
#' This function is designed for internal use and is not exported.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' w <- c(1, 1, 2, 2, 4)
#' weighted_quantiles(x, w, probs = c(0.25, 0.5, 0.75))
#'
#' @keywords internal
weighted_quantiles <- function(
    x,
    weights,
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95)
) {
  stopifnot(is.matrix(x), is.matrix(weights))
  stopifnot(all(dim(x) == dim(weights)))
  stopifnot(all(probs >= 0 & probs <= 1))
  
  n_rows <- nrow(x)
  n_probs <- length(probs)
  result <- matrix(NA_real_, nrow = n_rows, ncol = n_probs)
  colnames(result) <- paste0("q", probs * 100)
  
  # Collapse duplicate weights by averaging associated x values
  collapse_duplicates <- function(w, x) {
    uw <- unique(w)
    agg_x <- vapply(uw, function(wi) {
      mean(x[w == wi])
    }, numeric(1))
    list(w = uw, x = agg_x)
  }
  
  for (i in seq_len(n_rows)) {
    xi <- x[i, ]
    wi <- weights[i, ]
    
    # Remove NA values
    valid <- !(is.na(xi) | is.na(wi))
    xi <- xi[valid]
    wi <- wi[valid]
    
    if (length(xi) == 0 || sum(wi) == 0) next
    
    # Order by xi
    ord <- order(xi)
    xi <- xi[ord]
    wi <- wi[ord]
    
    # Normalize weights
    w_cum <- cumsum(wi) / sum(wi)
    
    # Collapse duplicated cumulative weights
    if (any(duplicated(w_cum))) {
      collapsed <- collapse_duplicates(w_cum, xi)
      w_cum <- collapsed$w
      xi <- collapsed$x
    }
    
    # Compute quantiles via interpolation
    qvals <- numeric(n_probs)
    for (j in seq_along(probs)) {
      p <- probs[j]
      if (p <= min(w_cum)) {
        qvals[j] <- xi[1]
      } else if (p >= max(w_cum)) {
        qvals[j] <- xi[length(xi)]
      } else {
        qvals[j] <- approx(w_cum, xi, xout = p)$y
      }
    }
    
    result[i, ] <- qvals
  }
  
  return(result)
}