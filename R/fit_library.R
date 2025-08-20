#' @title Build predictive/fucntional libraries
#' @description
#' This function builds a library of predictive functions based on memory-based learning
#' @usage
#' fit_library(Xr, Yr, sm, k, k_diss, k_range, ws = NULL, pls_c, 
#'             pc_selection = list("opc", min(dim(Xr), 40)),
#'             group, build_around, gh = TRUE, center = TRUE,
#'             scale = FALSE, diss_predictors = FALSE,
#'             metric = "rmse", return_best = TRUE,
#'             pls_max_iter = 1, pls_tol = 1e-6,
#'             documentation = character(), ...)
#' @param Xr a numeric matrix (or data.frame) of predictor variables of dimentions `n * p` corresponding to the reference data (observations in rows and variables in columns). 
#' @param Yr a numeric vector of length `n` containing the values of the response variable corresponding to the reference data. Missing values might be allowed (see details).
#' @param diss_method Either a character string indicating the dissimilarity
#'   metric to be used for selecting the nearest neighbors, or a precomputed
#'   dissimilarity matrix.
#'
#'   If a character is provided, it must be one of the following:
#'
#'   - `"pca"`: Mahalanobis distance based on principal component (PC) scores
#'     computed via singular value decomposition (SVD). See [ortho_diss()].
#'
#'   - `"pca.nipals"`: Same as `"pca"` but using the NIPALS algorithm. See
#'     [ortho_diss()].
#'
#'   - `"pls"`: Mahalanobis distance based on PLS scores. Requires `Yr`. See
#'     [ortho_diss()].
#'
#'   - `"mpls"`: Mahalanobis distance based on modified PLS scores (Shenk and
#'     Westerhaus, 1991; Westerhaus, 2014). Requires `Yr`. See [ortho_diss()].
#'
#'   - `"cor"`: Dissimilarity based on correlation between observations.
#'     See [cor_diss()].
#'
#'   - `"euclid"`: Euclidean distance. See [f_diss()].
#'
#'   - `"cosine"`: Cosine distance. See [f_diss()].
#'
#'   - `"sid"`: Spectral information divergence. See [sid()].
#'
#'   These options are passed to the [dissimilarity()] function, which performs
#'   the actual computation. Additional arguments to `dissimilarity()` can be
#'   passed through `...`.
#'
#'   If a matrix is provided, it must represent a precomputed dissimilarity
#'   structure with constraints depending on whether `anchor_indices` is used:
#'
#'   - If `anchor_indices = NULL` (default), `diss_method` must be a square
#'     matrix of size `n x n`, where `n = nrow(Xr)`. It must be symmetric
#'     with zeros on the diagonal, representing self-dissimilarity.
#'
#'   - If `anchor_indices` is specified (with length `m`), then
#'     `diss_method` must have dimensions `n x m`, where `n = nrow(Xr)` and
#'     `m = length(anchor_indices)`. The submatrix `diss_method[anchor_indices, ]`
#'     must be symmetric with zeros on the diagonal.
#'
#'   In all cases, the dissimilarity matrix must be consistent with the order
#'   of samples in `Xr`. That is, the rows of `diss_method` must correspond
#'   exactly to the rows of `Xr`, and when used, the elements of
#'   `anchor_indices` must refer to valid rows in both `Xr` and `diss_method`.
#' @param k a vector of integers specifying the sequence of k nearest neighbours to be tested. Either `k` or `k_diss` must be specified. This `vector` will be automatically sorted into ascending order. The values of vectors of class numeric, will be coerced to the next upper inetegers.
#' @param k_diss NOT YET IMPLEMENTED a numeric `vector` specifying the sequence of dissimilarity thresholds to be tested for the selection of the nearest neighbors around each observation. For a given observation, its neighbors are those that exhibit a dissimilarity value equal to or below a given threshold. These thresholds depend on the corresponding dissimilarity measure specified in `sm`. Either `k` or `k_diss` must be specified. 
#' @param k_range NOT YET IMPLEMENTED a vector of two integers specifying the minimum (first value) and the maximum (second value) number of neighbours allowed when the `k_diss` argument is used. 
#' @param ws an optional odd integer value which specifies a window size when the correlation dissimilarity is used (i.e `sm = "cor"`). If not specified, the standard correlation dissimilarity is computed (i.e no moving window). Default is `NULL`.
#' @param pls_c a vector of two integers indicating the minimum (first value) and maximum (second value) number of pls components to be used in the weighted average pls (wapls) regressions. 
#' @param pc_selection a list specifying the details of the method to be used for identifying the number of components to be retained for computing the dissimilarities between samples when the method passed to `sm` is one of the following options: `"pc"`, `"loc.pc"`, ` "pls"` or ` "loc.pls"`. This list must contain two elements in the following order: \itemize{
#'        \item{`method`:}{the method for selecting the number of components. Possible options are:  `"opc"` (optimized pc selection based on Ramirez-Lopez et al. (2013a, 2013b). See the \code{\link{orthoProjection}} function for more details;  `"cumvar"` (for selecting the number of principal components based on a given cumulative amount of explained variance); `"var"` (for selecting the number of principal components based on a given amount of explained variance); and  `"manual"` (for specifying manually the desired number of principal components)}
#'        \item{`value`:}{a numerical value that complements the selected method. If `"opc"` is chosen, it must be a value indicating the maximal number of principal components to be tested (see Ramirez-Lopez et al., 2013a, 2013b). If `"cumvar"` is chosen, it must be a value (larger than 0 and below 1) indicating the maximum amount of cumulative variance that the retained components should explain. If `"var"` is chosen, it must be a value (larger than 0 and below 1) indicating that components that explain (individually) a variance lower than this threshold must be excluded. If `"manual"` is chosen, it must be a value specifying the desired number of principal components to retain.
#'        }}
#' The default method for the `pc_selection` argument is `"opc"` and the maximum number of principal components to be tested is set to `min(dim(Xr), 40)`, where code{Xr} is the matrix of reference predictors in the \code{\link{mbl}} function.
#' Optionally, the `pc_selection` argument admits `"opc"` or `"cumvar"` or `"var"` or `"manual"` as a single character string. In such a case the default for `"value"` when either `"opc"` or `"manual"` are used is `min(dim(Xr), 40)`. When `"cumvar"` is used the default `"value"` is set to 0.99 and when `"var"` is used the default `"value"` is set to 0.01.
#' @param group (PUT THE DETAILS OF THE 1_NN VALIDATION IN DETAILS) an optional `factor` (or `vector` that can be coerced to a `factor` by `as.factor`) that assigns to each observation in `Xr` a group/class label (e.g. groups can be given by spectra collected from the same batch of measurements, from the same sample, from samples with very similar origin, etc). This is taken into account for internal validation of pls models (and factor optimization) to avoid pseudo-replication. For validation, the model of a given target observation is first fitted based on its neighbors and excluding that observation (and all the samples belonging to its group), then the model is used to predict the response variable of the target observation. See details.
#' @param anchor_indices An optional integer vector of length `m` (`m < n`) specifying the row indices of the reference samples around which local models should be constructed. If `NULL` (default), local models will be built for all `n` samples in the reference set. See Details.
#' @param gh a logical indicating whether or not to compute and return the Mahalanobis distance (in the pls space) between each element in `Xr` and the center of `Xr`.
#' @param center a logical indicating whether or not the predictor variables must be centered at each local segment (before regression).
#' @param scale a logical indicating whether or not the predictor variables must be scaled to unit variance at each local segment (before regression).
#' @param diss_predictors a logical indicating if the local (square symmetric) dissimilarity matrix corresponding the selected neighbourhood shall be used as source of additional predictors (i.e the columns of this local matrix are treated as predictor variables). In some cases this may result in an improvement  of the prediction performance (Ramirez-Lopez et al., 2013a). 
#' @param metric a character value indicating what model perfoamance statistics shall be used to select the best set of parameters (i.e. number of neighbors or threshold distances and pls factors)  to fit the final model. Options are `"rmse"` (default) or `"r2"` (coefficient of determination).
#' @param return_best a logical indicating if the final library of functions using the optimal parameters found shall be returned.
#' @param pls_max_iter (BETTER DESCRIPTION REQUIRED) maximum number of iterations for the partial least squares methods.
#' @param pls_tol (BETTER DESCRIPTION REQUIRED) for convergence in the partial orthogonal scores partial least squares regressions using the nipals algorithm. Default is 1e-6
#' @param chunk_size Integer. Specifies the maximum number of samples
#'   (rows or columns, depending on context) to include in each chunk.
#'   Each chunk is processed independently, which can be beneficial
#'   for parallel computation. For example, each chunk can be dispatched
#'   to a separate processor or core, where the samples within the chunk
#'   are processed sequentially. Defaults to `1` (i.e., one sample per
#'   chunk). The value of `chunk_size` cannot exceed the total number of
#'   available samples for the given context (e.g., `nrow(Xr)` or 
#'   `length(anchor_indices)` when provided).
#' @param documentation (BETTER DESCRIPTION REQUIRED) an optional character string for documentating the call to this function.
#' @param ...  arguments passed to the \code{\link{dissimilarity}} function.
#' @details
#' 
#' ** Anchor indices **  
#' 
#' By default, local models are constructed around all `n` samples in the
#' reference set. Alternatively, the user may specify a subset of `m` samples
#' (`m < n`) around which the library of local models should be built. This is
#' done by supplying their row indices to the `anchor_indices` argument.
#'
#' In this configuration, each local model is still constructed using the
#' nearest neighbors selected from the full set of `n` reference samples, but
#' only for the `m` anchor samples. This approach is especially useful when the
#' reference set is very large, as it reduces the total number of models while
#' retaining representativeness and predictive performance—particularly if the
#' anchor samples are chosen to reflect the diversity of the dataset.
#'
#' This strategy can lead to substantial computational savings and a more
#' compact library of models. **Note** that when distance calculations rely on
#' the response variable (e.g., GH distance, PLS distance, or optimized PC
#' distance), the response values (`Yr`) of the anchor samples themselves are
#' not used during distance computation. This is done to improve computational
#' efficiency: dissimilarities are computed between each anchor sample and all
#' reference samples (excluding the anchor itself), rather than recalculating
#' the full `n × n` distance matrix. This significantly reduces the number of
#' response-based distance computations when the number of anchors is small.
#' Also, **note** that the response values (Yr) of the anchor observations are 
#' always used during the fitting of  their corresponding local models.
#' For efficiency, the number of anchor indices must be less than 90% of the
#' total number of rows in `Xr`. If this threshold is exceeded, consider
#' building a model for each sample in `Xr` by setting `anchor_indices = NULL`.
#'
#' ** Missing values in `Yr` **  
#' 
#' Missing values in `Yr` are allowed. If they exceed 25% of the total number
#' of observations, a warning will be issued. The local model of an observation
#' with a missing `Yr` value is fitted only using its nearest neighbors (i.e.,
#' without including the observation itself).
#' @references 
#' Rajalahti, T., Arneberg, R., Berven, F. S., Myhr, K. M., Ulvik, R. J., & 
#' Kvalheim, O. M. (2009). Biomarker discovery in mass spectral profiles by 
#' means of selectivity ratio plot. Chemometrics and Intelligent Laboratory 
#' Systems, 95(1), 35-48.
#'
#' @return DESCRIOTION REQUIRED
#' @author Leonardo Ramirez-Lopez
#' @examples
#' \dontrun{
#' #' ## GOOD EXAMPLES ARE REQUIRED
#' }
#' @export fit_library

fit_library <- function(
    Xr,
    Yr,
    diss_method,
    k,
    k_diss,
    k_range,
    ws = NULL,
    pls_c,
    pc_selection = list(method = "var", value = 0.01),
    group,
    anchor_indices = NULL,
    gh = TRUE,
    center = TRUE,
    scale = FALSE,
    diss_predictors = FALSE,
    metric = "rmse",
    return_best = TRUE,
    pls_max_iter = 1,
    pls_tol = 1e-6,
    nearest_neighbor_val= TRUE,
    optimize_ncomp_range = FALSE,
    chunk_size = 1,
    documentation = character(),
    ...
) {
  
  call.f <-(match.call())    
  
  if (is.null(colnames(Xr))) {
    stop("column names are mandatory for Xr")
  }
  
  if (missing(k)) {
    stop("'k' is missing... at the mooment 'k_diss' is not active")
  }
  
  
  if (missing(group)) {
    group <- as.factor(paste("g", 1:nrow(Xr), sep = ""))
  }
  ## Compute the dissimilarity matrix for Xr
  
  if (class(Yr) != "matrix") {
    Yr <- matrix(Yr, nrow = nrow(Xr))
  }
  
  
  if (missing(k) && missing(k_diss)) {
    stop("Either k or k_diss must be specified")
  }
  
  if (!missing(k) && !missing(k_diss)) {
    stop("Only one of k or k_diss can be specified")
  }
  
  if (!missing(k)) {
    if (!is.numeric(k)) {
      stop("k must be a numeric vector")
    }
    k <- sort(unique(ceiling(k)))
  }
  
  if (!missing(k_diss)) {
    dtc <- k_diss
    if (missing(k_range)) {
      stop("When using k_diss, 'k_range' must be specified")
    }
    
    if (!is.numeric(k_range) || length(k_range) != 2 || any(is.na(k_range))) {
      stop("'k_range' must be a numeric vector of length 2")
    }
    
    k_range <- ceiling(k_range)
    
    if (k_range[1] > k_range[2]) {
      stop("In 'k_range', the first value must be <= the second value")
    }
    
    k.min <- as.integer(k_range[1])
    k.max <- as.integer(k_range[2])
    if (k.min < 10)
      stop("Minimum number of nearest neighbours allowed is 10")
    if (k.max > nrow(Xr))
      stop("Maximum number of nearest neighbours cannot exceed the number of reference observations")
  } else{
    dtc <- NULL
  }
  
  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size < 1) {
    stop("chunk_size must be a single numeric value larger than 0")
  }
  
  min_size <- if (!is.null(anchor_indices)) {
    min(nrow(Xr), length(anchor_indices))
  } else {
    nrow(Xr)
  }
  
  if (chunk_size > min_size) {
    stop("chunk_size cannot exceed ", min_size, " in this context")
  }
  
  pc_selection <- check_pc_arguments(
    n_rows_x = nrow(Xr), 
    n_cols_x = ncol(Xr),  
    pc_selection = pc_selection
  )$pc_selection_checked
  
  if (!is.null(anchor_indices)) {
    if (!is.numeric(anchor_indices) || any(anchor_indices < 1) || any(anchor_indices > nrow(Xr))) {
      stop("anchor_indices must be a numeric vector with valid row indices of Xr")
    }
    if (length(anchor_indices) > floor(nrow(Xr) * 0.90)) {
      stop(
        paste0(
          "Too many anchor indices specified: must be less than 90% of the rows ",  
          "in Xr. Consider building a model for each sample in Xr by ", 
          "setting 'anchor_indices = NULL'."
        )
      )
    }
  }
  
  if (is.matrix(diss_method)) {
    diss_method_type <- "Precomputed dissimilarity matrix"
    if (nrow(diss_method) != nrow(Xr)) {
      stop(
        "If 'diss_method' is a matrix, it must have the same number of rows ",
        "as 'Xr'."
      )
    }
    
    if (!is.null(anchor_indices)) {
      if (ncol(diss_method) != length(anchor_indices)) {
        stop(
          "If 'diss_method' is a matrix and 'anchor_indices' is provided, ",
          "it must have the same number of columns as the length of ",
          "'anchor_indices'."
        )
      }
      
      if (any(abs(diag(diss_method[anchor_indices, ])) > 1e-8)) {
        stop(
          "'diss_method[anchor_indices, ]' must be symmetric with zeros on ", 
          "the diagonal. Ensure the dissimilarity matrix is correctly defined ", 
          "for the selected anchor samples."
        )
      }
      
    } else {
      if (ncol(diss_method) != nrow(Xr)) {
        stop(
          "If 'diss_method' is a matrix and 'anchor_indices' is NULL, ",
          "it must be square with dimensions equal to 'nrow(Xr)'."
        )
      }
      
      if (any(abs(diag(diss_method)) > 1e-8)) {
        stop(
          "The diagonal of 'diss_method' must contain only zeros ",
          "(i.e., self-dissimilarities must be zero)."
        )
      }
    }
    dsm <- list()
    dsm$dissimilarity <- diss_method
    rm(diss_method)
    gc()
    sml <- list(diss_method = "Precomputed dissimilarity matrix")
    dsm <- append(sml, dsm)
    if (!is.null(anchor_indices)) {
      Yr_anchor <- Yr[anchor_indices, , drop = FALSE]
    } else {
      Yr_anchor <- Yr
    }
  } else if (is.character(diss_method)) {
    diss_method_type <- match.arg(
      diss_method,
      c("pca", "pca.nipals", "pls", "mpls", "cor", "euclid", "cosine", "sid")
    )
    
    cat("Computing dissimilarities... \n")
    
    if (!is.null(anchor_indices)) {
      Yr_anchor <- Yr[anchor_indices, , drop = FALSE]
      if (sum(!is.na(Yr_anchor)) < 3) {
        stop("At least 3 non-missing values in Yr are required for the anchor indices")
      }
      dsm <- dissimilarity(
        Xr = Xr[-anchor_indices, , drop = FALSE],
        Xu = Xr[anchor_indices, , drop = FALSE],
        diss_method = diss_method_type,
        Yr = Yr[-anchor_indices, , drop = FALSE],
        center = center,
        scale = scale,
        gh = gh,
        pc_selection = pc_selection,
        return_projection = TRUE,
        ws = ws
      )
      if (gh) {
        gh_c <- rep(NA, nrow(Xr))
        gh_c[-anchor_indices] <- dsm$gh$gh_Xr
        gh_c[anchor_indices] <- dsm$gh$gh_Xu
        dsm$gh$gh_Xr <- gh_c
        dsm$gh <- dsm$gh[!names(dsm$gh) %in% "gh_Xu"]
      } else {
        dsm$gh <- NULL
      }
      
      new_mat <- matrix(NA_real_, nrow(Xr), ncol(dsm$dissimilarity))
      new_mat[-anchor_indices, ] <- dsm$dissimilarity
      
      # Assign
      dsm$dissimilarity <- new_mat
      rm(new_mat)
      gc() 
      
      if (diss_method_type %in% c("pca", "pca.nipals", "pls", "mpls")) {
        z_anchor <- dsm$projection$scores[-(1:(nrow(Xr) - length(anchor_indices))), ] 
        z_anchor <- scale(
          z_anchor, 
          center = FALSE, 
          scale = dsm$projection$scores_sd
        )
        z_diss <- fast_self_euclid(z_anchor)
        dsm$dissimilarity[anchor_indices, ] <- z_diss
        rm(z_diss, z_anchor)
        gc() 
        real_order <- order(c((1:nrow(Xr))[-anchor_indices], anchor_indices))
        dsm$projection$scores <- dsm$projection$scores[real_order, ]
        rownames(dsm$projection$scores) <- paste0("Xr_", 1:nrow(Xr))
        
        dsm$projection$scores[real_order, ][anchor_indices, ] |> rownames()
        
      } else {
        z_anchor <- Xr[anchor_indices, , drop = FALSE]
        if (center || scale) {
          z_anchor <- scale(
            z_anchor,
            center = if (center) get_column_means(Xr) else FALSE,
            scale  = if (scale)  get_column_sds(Xr)  else FALSE
          )
        }
        dsm_anchor <- dissimilarity(
          Xr = Xr[anchor_indices, , drop = FALSE],
          diss_method = diss_method_type,
          Yr = Yr[anchor_indices, , drop = FALSE],
          center = FALSE,
          scale = FALSE,
          gh = FALSE,
          pc_selection = NULL,
          return_projection = FALSE,
          ws = ws
        )
        
        dsm$dissimilarity[anchor_indices, ] <- dsm_anchor$dissimilarity
        rm(dsm_anchor, z_anchor)
        gc() 
      }
    } else {
      Yr_anchor <- Yr
      dsm <- dissimilarity(
        Xr = Xr,
        diss_method = diss_method_type,
        Yr = Yr,
        center = center,
        scale = scale,
        gh = gh,
        pc_selection = pc_selection,
        return_projection = TRUE,
        ws = ws
      )
    }
    
    dsm <- dsm[!names(dsm) %in% "documentation"]
    sml <- list(diss_method = diss_method_type)
    dsm <- append(sml, dsm)
  } else {
    stop("diss_method must be a character string or a precomputed matrix")
  }
  
  max_k <- ifelse(missing(k), max(k_range), max(k))
  if (anyNA(Yr)) {
    if (sum(!is.na(Yr)) < max_k) {
      stop(
        paste0(
          "The maximum number of neighbors selected (", max_k, 
          ") exceeds the number of non-missing observations in 'Yr' (", 
          sum(!is.na(Yr)), "). Please reduce the number of neighbors."
        )
      )
    }
  }
  
  ## find the indices of the nearest neighbors
  kidxmat <- top_k_order(
    dsm$dissimilarity, k = max_k, skip = which(is.na(Yr))
  )
  
  kdissmat <- extract_by_index(dsm$dissimilarity, kidxmat)
  
  
  ## for each sample in Xu show what of its
  ## nearest neighbors samples belong to its group
  if (is.null(group)) {
    kidxgrop <- matrix(TRUE, nrow = nrow(kidxmat), ncol = ncol(kidxmat))
    kidxgrop[1, ] <- FALSE
  } else {
    kidxgrop <- not_in_same_group(kidxmat, group = group)
  }
  
  probs <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.00)
  dnms <- list(
    1:ncol(kidxgrop), 
    paste0(probs * 100, "%")
  )
  nnstats <- vector("list", length(k))
  for (ii in seq_along(k)) {
    nnstats[[ii]] <- compute_nn_quantiles(
      kidxmat = kidxmat,
      kidxgrop = kidxgrop,
      Yr = Yr,
      k = k[ii],
      probs = probs
    ) 
    dimnames(nnstats[[ii]]) <- dnms
  }
  names(nnstats) <- paste0("k.", k)
  
  
  ## RMSE/(q"75%" - q"25%")
  ## Similar to the ratio of performance to 
  ## inter-quartile distance (RPIQ)
  # Bellon-Maurel, V., Fernandez-Ahumada, E., 
  # Palagos, B., Roger, J.M., McBratney, A., 2010. 
  # Critical review of chemometric indicators commonly 
  # used for assessing the quality of the prediction of 
  # soil attributes by NIR spectroscopy. TrAC Trends 
  # in Analytical Chemistry, 29(9), pp.1073-1081.
  itq <- vapply(nnstats, function(x) {
    x[, "75%"] - x[, "25%"]
  }, numeric(nrow(nnstats[[1]])))
  
  if (diss_predictors) {
    dssm <- dsm$dissimilarity
  } else {
    dssm <- NULL
    npredictors <- ncol(Xr)
  }
  
  addit <- ifelse(return_best, 1, 0)
  
  pb <- txtProgressBar(min = 0, max = length(k) + addit, char = "-")
  
  minF <- min(pls_c)
  maxF <- max(pls_c)
  
  if (optimize_ncomp_range) {
    # what this grid does is to describe the possible combinations of min factor 
    # and max factors, this for optimizing the min and max pls factors 
    sgrid <- expand.grid(minpls = minF:maxF, maxpls = minF:maxF)
    sgrid <- sgrid[sgrid$minpls <= sgrid$maxpls, ]
    row.names(sgrid) <- 1:nrow(sgrid)
    
    # here the combinations of min and max pls factors is translated into a binary 
    # grid to switch on and off some factoes
    emgrid <- t(
      sapply(
        1:nrow(sgrid),
        FUN = function(q, x, wv) {
          wv[x[q, 1]:x[q, 2]] <- 1
          wv
        },
        x = sgrid,
        wv = rep(0, maxF)
      )
    )
    emgrid <- emgrid[, minF:maxF]
  } else {
    emgrid <- matrix(1, nrow= 1, length(minF:maxF))
    sgrid <- data.frame(
      minpls = minF,
      maxpls = maxF
    )
  }
  
  
  cat("Fitting models... \n")
  ## perform the nearest neighbor predictions  
  
  if (nearest_neighbor_val) {
    
    nnpreds <- sapply(
      1:length(k), 
      FUN = .get_all_fits,
      Xr = Xr, 
      Yr = Yr, 
      k = k,
      min_component = minF, 
      max_component = maxF, 
      emgrid = emgrid,
      scale = scale, 
      maxiter = pls_max_iter, 
      tol = pls_tol, 
      regression = TRUE, 
      pc_selection = pc_selection,
      kidxmat = kidxmat,
      kidxgrop = kidxgrop,
      dissimilarity_mat = dssm,
      pb = pb,
      chunk_size = chunk_size
    )
    
    ## Organize the results (in nnpreds)
    pparam <- matrix(NA, nrow(emgrid), 4)
    sstats <- function(y, yhat, itqk, pparam) {
      me <- sweep(-yhat, MARGIN = 1, STATS = y, FUN = "+", check.margin = FALSE)
      pparam[, 1] <- cor(y, yhat, use = "complete.obs")^2
      pparam[, 2] <- colMeans(me^2, na.rm = TRUE)^0.5
      pparam[, 3] <- colMeans(me, na.rm = TRUE)
      pparam[, 4] <- colMeans(
        sweep(me, 
              MARGIN = 1, 
              STATS = itqk, 
              FUN = "/", 
              check.margin = FALSE)^2, 
        na.rm = TRUE
      )^0.5
      return(pparam)
    }
    
    isubset3Row <- function(x1, x2, x3, x4) {
      sargs <- names(match.call())[-1]
      nextEl <- function(..ii..) {
        sapply(sargs, FUN = function(x) get(x)[,..ii..], simplify = FALSE)
      }
      obj <- list(nextElem = nextEl)
    }
    
    itr <- isubset3Row(x1 = nnpreds, x2 = itq)
    
    kpredstats <- function(
    ..k.., 
    itr,
    pparam,
    y
    ) {
      ne <- itr$nextElem(..k..)
      statsresults <- sstats(
        y = y, 
        yhat = t(matrix(ne$x1, nrow(pparam))), 
        itqk = ne$x2, 
        pparam = pparam
      )  
      statsresults
    }
    
    predperformance <- lapply(
      1:ncol(nnpreds), 
      FUN = kpredstats,
      itr = itr,
      pparam = pparam,
      y = Yr_anchor
    )
    
    predperformance <- data.frame(do.call("rbind", predperformance))
    colnames(predperformance) <- c("r2", "rmse", "me", "st.rmse")
    
    predperformance <- data.frame(
      minpls = rep(sgrid$minpls, times = length(k)),
      maxpls = rep(sgrid$maxpls, times = length(k)),
      k = rep(k, each = nrow(pparam)),
      predperformance
    )
    
    
    # ## store results
    # kresults <- data.frame(k = k, 
    #                        r2 = cor(nnpreds , Yr[,], use = "complete.obs")^2,
    #                        rmse = colMeans(me^2, na.rm = TRUE)^0.5,
    #                        rmse.st = colMeans((me/itq)^2, na.rm = TRUE)^0.5,
    #                        check.names = FALSE)
    # plot(kresults$k, kresults$rmse)
    
    # find optinmal parameters
    
    if (metric == c("rmse")) {
      bestp <- predperformance[which.min(predperformance[, metric]), ][1, ]
    }
    
    if (metric == c("r2")) {
      bestp <- predperformance[which.max(predperformance[, metric]), ][1, ]
    }
    
    optimalk <- bestp$k
    optimalminpls <- bestp$minpls
    optimalmaxpls <- bestp$maxpls
    
    
    ## Extract the vector of predictions corresponding to the best predictions
    ## (optimal k, optimal pls range)
    plsitemn <- which(
      sgrid$minpls == optimalminpls & sgrid$maxpls == optimalmaxpls
    )
    plsitemn <- seq(plsitemn, by = nrow(pparam), length.out = length(Yr_anchor))
    bestpreds <- itr$nextElem(which(k == optimalk))$x1[plsitemn]
    bestpredsresiduals <- Yr_anchor - bestpreds   
    
    # ithbarrio <- ith_subsets(x = Xr,
    #                          y = Yr,
    #                          kindx = kidxmat[1:optimalk,],
    #                          D = dssm)
  } else {
    optimalk <- max(k)
    optimalminpls <- min(pls_c)
    optimalmaxpls <- max(pls_c)
    bestpredsresiduals <- bestp <- predperformance <- NULL
    setTxtProgressBar(pb, 1)
  }

  if (return_best) {
    ## build the pls librarby
    # plslib <- sapply(1:nrow(Xr), 
    #                  FUN = final_fits,
    #                  ithsubset = ithbarrio,
    #                  minF = optimalminpls, 
    #                  maxF = optimalmaxpls, 
    #                  scale = group, 
    #                  maxiter = pls_max_iter, 
    #                  tol = pls_tol, 
    #                  pc_selection = pc_selection)
    
    if (diss_predictors) {
      n_var <- 1 + (5 * (optimalk + ncol(Xr)))
    } else {
      n_var <- 1 + (5 * ncol(Xr))
    }
    plslib_template <- matrix(
      NA, 
      nrow = n_var, 
      ncol = chunk_size
    )
    n_iter <- ceiling(ncol(dsm$dissimilarity) / chunk_size)
    plslib <- foreach(
      i = 1:n_iter, 
      .export = c(
        "ith_pred_subsets",
        "ith_subsets_list",
        "ith_subsets_by_group",
        ".get_all_fits",
        "ith_local_fit", 
        "final_fits_cpp"
      ),
      ithbarrio = ith_subsets_list(
        x = Xr, y = Yr, kindx = kidxmat[1:optimalk, ], D = dssm, chunk_size = chunk_size
      )
    ) %dopar% { 
      iplslib <- plslib_template
      for (j in 1:length(ithbarrio)) {
        ij_pls <- final_fits_cpp(
          X = ithbarrio[[j]]$x,
          Y = ithbarrio[[j]]$y,
          new_x = ithbarrio[[j]]$x[1, , drop = FALSE],
          min_component = optimalminpls, 
          max_component = optimalmaxpls, 
          scale = scale, 
          maxiter = pls_max_iter, 
          tol = pls_tol
        )
        
        # ij_pls <- final_fits(
        #   ilocalsubset = ithbarrio,
        #   min_component = optimalminpls, 
        #   max_component = optimalmaxpls, 
        #   scale = scale, 
        #   maxiter = pls_max_iter, 
        #   tol = pls_tol
        # )
        iplslib[, j] <- unlist(ij_pls, recursive = FALSE, use.names = FALSE)
      }
      if (j < ncol(iplslib))  {
        iplslib <- iplslib[, 1:j]
      }
      iplslib
    }
    
    plslib <- do.call("cbind", plslib)
    
    setTxtProgressBar(pb, length(k) + 1)
    
    plslib <- t(plslib)
    
    namesk <- NULL
    if (diss_predictors) {
      npredictors <- optimalk + ncol(Xr)
      namesk <- paste("k", 1:optimalk, sep = "")
    }
    if (!diss_predictors) {
      npredictors <- ncol(Xr)
    }
    
    xscale <- plslib[ , -c(1:((4 * npredictors) + 1))]
    plslib <- plslib[ , c(1:((4 * npredictors) + 1))]
    
    xcenter <- plslib[ , -c(1:((3 * npredictors) + 1))]    
    plslib <- plslib[ , c(1:((3 * npredictors) + 1))]
    
    plsvips <- plslib[ ,-c(1:(npredictors + 1), (ncol(plslib) - npredictors + 1):ncol(plslib))]
    plssratios <- plslib[ ,c((ncol(plslib) - npredictors + 1):ncol(plslib))]
    plslib <- plslib[ , c(1:(npredictors + 1))]
    
    colnames(plslib) <- c("b0", namesk, colnames(Xr))
    colnames(xcenter) <- colnames(xscale) <- c(namesk, colnames(Xr))
    colnames(plssratios) <- colnames(plsvips) <- c(namesk, colnames(Xr))
    
    if (diss_predictors) {
      bs <- list(
        B0 = plslib[,1],
        B = plslib[,colnames(plslib) %in% colnames(Xr)],
        Bk = plslib[,colnames(plslib) %in% namesk]
      )
    } else {
      bs <- list(
        B0 = plslib[,1],
        B = plslib[,colnames(plslib) %in% colnames(Xr)]
      )
    }
    
    if (center) {
      center <- get_column_means(Xr)
    } else{
      center <- rep(0, ncol(Xr))
    }
    
    if (scale) {
      gscale <- get_column_sds(Xr)
    } else {
      gscale <- rep(1, ncol(Xr))
    }
    
    iscale <- list(
      centre = center,
      scale = gscale,
      local.x.center = xcenter[ , colnames(xscale) %in% colnames(Xr)],
      local.x.scale = xscale[ , colnames(xscale) %in% colnames(Xr)]
    )
    if (diss_predictors) {
      iscale$local.diss.scale <- xscale[ ,colnames(xscale) %in% namesk]
    }
    
    if (diss_method_type %in% c("pca", "pls")) {
      fresults <- list(
        dissimilatity = dsm[!names(dsm) %in% "gh"],
        gh = dsm$gh, 
        results = predperformance,
        best = bestp,
        sel.param = list(
          optimalk = optimalk,
          optimal.factors = c(min.pls = bestp$minpls, max.pls = bestp$maxpls)
        ),
        residuals = bestpredsresiduals,
        functionlibrary = bs,
        functionvips = plsvips,
        functionselectivityrs = plssratios,
        scale =  iscale,
        yu.nnstats = nnstats,
        documentation = documentation
      )
    } else {
      fresults <- list(
        dissimilatity = dsm[!names(dsm) %in% "gh"],
        gh = dsm$gh, 
        results = predperformance,
        best = bestp,
        sel.param = list(
          optimalk = optimalk,
          optimal.factors = c(min.pls = bestp$minpls, max.pls = bestp$maxpls)
        ),
        residuals = bestpredsresiduals,
        functionlibrary = bs,
        functionvips = plsvips,
        functionselectivityrs = plssratios,
        scale = iscale,
        # Xr = Xr,
        yu.nnstats = nnstats,
        documentation = documentation
      )
    }
    
    if (!is.null(rownames(Xr))) {
      namesrows <- rownames(Xr)
    } else{
      namesrows <- 1:nrow(Xr)
    }
    
    if (!is.null(anchor_indices)) {
      namesrows <- namesrows[anchor_indices]
    }
    
    fresults$yu.nnstats <- lapply(
      fresults$yu.nnstats, 
      FUN = function(x, nms) {
        rownames(x) <- nms
        x
      }, 
      nms = namesrows
    )
    
    rownames(fresults$functionvips) <- 
      rownames(fresults$functionselectivityrs) <- 
      rownames(fresults$functionlibrary$B) <- 
      names(fresults$functionlibrary$B0) <- namesrows
    if (!is.null(fresults$functionlibrary$Bk))
      rownames(fresults$functionlibrary$Bk) <- namesrows
    class(fresults) <- c("funlib","list")
  } else {
    
    fresults <- list(
      dissimilatity = dsm[!names(dsm) %in% "gh"],
      gh = dsm$gh, 
      results = predperformance,
      best = bestp,
      residuals = bestpredsresiduals,
      yu.nnstats = nnstats,
      documentation = documentation
    )
    
    if (!is.null(rownames(Xr))) {
      namesrows <- rownames(Xr)
    } else{
      namesrows <- 1:nrow(Xr)
    }
    
    if (!is.null(anchor_indices)) {
      namesrows <- namesrows[anchor_indices]
    }
    
    fresults$yu.nnstats <- lapply(
      fresults$yu.nnstats, 
      FUN = function(x, nms) {
        rownames(x) <- nms
        x
      }, 
      nms = namesrows
    )
    class(fresults) <- c("validationfunlib","list")
  }
  fresults$anchor_indices <- anchor_indices
  
  if (nearest_neighbor_val) {
    rownames(fresults$residuals) <- namesrows
  }
  
  attr(fresults, "call") <- call.f
  return(fresults)
}



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
    residual_cutoff = NULL,
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
  )[1:object$sel.param$optimalk, ]
  

  ## this might become an argument to cancel models with high residuals (rd > xx)
  if (!is.null(object$residuals) & !is.null(residual_cutoff)) {
    abs_res <- abs(object$residuals)
    # res <- abs(scale(res, center = TRUE, scale = TRUE))
    plot(abs_res)
    abs_res[abs_res < residual_cutoff] <- 0
    abs_res[abs_res >= residual_cutoff] <- 1
    high_residual_model <- as.logical(abs_res)
    high_residual_model[is.na(high_residual_model)] <- TRUE
    
    xudss <- sapply(
      1:ncol(xunn),
      FUN = function(diss, nn, cs, index){
        d1 <- diss[nn[, index], index]
        d2 <- cs[nn[, index]]
        # assign the max dissimilarity to make sure the model is not selected
        d1[d2] <- max(d1) 
        d1
      },
      diss = dsmxu$dissimilarity, 
      nn = xunn,
      cs = high_residual_model
    )[1:object$sel.param$optimalk, ]
  } else {
    xudss <- apply(
      dsmxu$dissimilarity, 
      MARGIN = 2, 
      FUN = sort
    )[1:object$sel.param$optimalk, ]
    
    high_residual_model <- rep(FALSE, ncol(xunn))
  }
  ## this deactivates model cancelling (for the moment)
  # res[] <- FALSE 
  
  dweights <- sweep(
    xudss, 
    MARGIN = 2, 
    STATS = get_column_maxs(xudss), 
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
  
  # Transpose weights 
  dweights <- t(dweights)
  
  weighted_yu_preds <- yupreds * dweights
  
  # Weighted mean per row
  wpredictions <- rowSums(weighted_yu_preds)
  
  # Centered deviations (from original, unweighted predictions)
  centered_dev <- sweep(yupreds, 1, wpredictions, FUN = "-")^2
  
  # Apply weights to squared deviations
  wvar <- rowSums(centered_dev * dweights)
  
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

