#' @title iterator for nearest neighbor subsets
#' @param x a reference matrix
#' @param xu a second matrix
#' @param y a matrix of side information
#' @param kindx a matrix of nearest neighbor indices
#' @param na_rm logical indicating whether NAs must be removed.
#' @description internal
#' @keywords internal
ith_subsets_ortho_diss <- function(x,
                                   xu = NULL,
                                   y,
                                   kindx,
                                   na_rm = FALSE) {
  it_kindx <- iter(kindx, by = "column")
  it_xu <- iter(xu, by = "row")

  nextEl <- function() {
    knns <- nextElem(it_kindx)
    if (na_rm) {
      knns <- knns[!is.na(knns)]
    }
    lsubst <- list(x = x[knns, ])
    if (!is.null(y)) {
      lsubst$y <- y[knns, , drop = FALSE]
    }
    if (!is.null(xu)) {
      lsubst$xu <- nextElem(it_xu)
    }
    return(lsubst)
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}

#' @title local ortho dissimilarity matrices initialized by a global
#' dissimilarity matrix
#' @param k_index_matrix a matrix of nearest neighnbor indices
#' @param Xr argument passed to ortho_projection
#' @param Yr argument passed to ortho_projection
#' @param Xu argument passed to ortho_projection
#' @param diss_method argument passed to ortho_projection
#' @param pc_selection argument passed to ortho_projection
#' @param center argument passed to ortho_projection
#' @param scale argument passed to ortho_projection
#' @description internal
#' @keywords internal
local_ortho_diss <- function(k_index_matrix, Xr, Yr, Xu,
                             diss_method, pc_selection, center, scale,
                             allow_parallel, ...) {
  "%mydo%" <- get("%do%")
  if (allow_parallel & getDoParRegistered()) {
    "%mydo%" <- get("%dopar%")
  }

  if (!is.null(Xu)) {
    n_xu <- nrow(Xu)
  } else {
    n_xu <- 0
  }

  n_xr <- nrow(Xr)

  if (diss_method %in% c("pca", "pca.nipals")) {
    use_distance_method <- "euclid"
  } else {
    use_distance_method <- "mahalanobis"
  }

  if (n_xr + n_xu == ncol(k_index_matrix)) {
    ## squared matrix for all possible pair-wise dissimilarities
    distnc <- matrix(NA, ncol(k_index_matrix), ncol(k_index_matrix))
    diag(distnc) <- 0
    ith_subsets <- ith_subsets_ortho_diss(
      x = Xr, y = Yr,
      kindx = k_index_matrix
    )
    n_rows_d <- ncol(k_index_matrix)
    # neighbors (excluding target sample )
    n_neighbors <- nrow(k_index_matrix) - 1
    has_target_indices <- TRUE
  } else {
    ## matrix for dissimilarities between Xr and Xu
    distnc <- matrix(NA, n_xr, ncol(k_index_matrix))
    ith_subsets <- ith_subsets_ortho_diss(
      x = Xr, xu = Xu, y = Yr,
      kindx = k_index_matrix,
      na_rm = TRUE
    )
    n_rows_d <- n_xr
    n_neighbors <- nrow(k_index_matrix)
    has_target_indices <- FALSE
  }


  dist_vec <- rep(NA, n_neighbors + 2)

  ith_local_subset <- NULL
  ## for each computations
  i <- NULL
  local_d <- foreach(
    i = 1:ncol(k_index_matrix),
    ith_local_subset = ith_subsets,
    .combine = cbind,
    .inorder = FALSE,
    .export = c(
      "ortho_projection",
      "pc_projection",
      "pls_projection",
      "sim_eval",
      "dist_vec"
    ),
    .packages = c("resemble")
  ) %mydo% { # Parallel
    ith_projection <- ortho_projection(
      Xr = ith_local_subset$x,
      Yr = ith_local_subset$y,
      Xu = ith_local_subset$xu,
      method = diss_method,
      pc_selection = pc_selection,
      center = center,
      scale = scale,
      ...
    )

    if (diss_method %in% c("pca", "pca.nipals")) {
      scores <- sweep(ith_projection$scores,
        MARGIN = 2,
        STATS = ith_projection$scores_sd,
        FUN = "/"
      )
    } else {
      scores <- ith_projection$scores
    }

    ith_dist <- dist_vec
    ith_dist[3:length(ith_dist)] <- f_diss(
      Xr = scores[1, , drop = FALSE],
      Xu = scores[-1, , drop = FALSE],
      diss_method = use_distance_method,
      center = FALSE,
      scale = FALSE
    )
    ith_dist[1] <- i
    ith_dist[2] <- ith_projection$n_components
    ith_dist
  }
  ## end of foreach computations


  local_d <- local_d[-1, order(local_d[1, ])]
  local_n_components <- local_d[1, ]
  local_d <- local_d[-1, ]
  add_row_indices <- ((0:(ncol(k_index_matrix) - 1)) * n_rows_d)
  # necessary to remove the xi-to-xi from the indices of nearest neighnbors
  if (has_target_indices) {
    k_index_matrix <- k_index_matrix[-1, ]
  }
  vec_positions <- sweep(k_index_matrix,
    MARGIN = 2,
    STATS = add_row_indices, FUN = "+"
  )
  distnc[vec_positions] <- as.vector(local_d)
  return(list(
    dissimilarity_m = distnc,
    local_n_components = data.table(local_n_components)
  ))
}



#' @title format internal messages
#' @param xr_xu_names the names of Xr and Xu
#' @description internal
#' @keywords internal
format_xr_xu_indices <- function(xr_xu_names) {
  xu_insufficient <- gsub("Xu_", "", xr_xu_names[grep("^Xu_", xr_xu_names)])
  xr_insufficient <- gsub("Xr_", "", xr_xu_names[grep("^Xr_", xr_xu_names)])


  xu_mss <- paste0(c(
    ifelse(length(xu_insufficient) > 0,
      "\nXu: ", ""
    ),
    paste(xu_insufficient, collapse = ", ")
  ),
  collapse = " "
  )
  xr_mss <- paste0(c(
    ifelse(length(xr_insufficient) > 0,
      "\nXr: ", ""
    ),
    paste(xr_insufficient, collapse = ", ")
  ),
  collapse = " "
  )

  return(list(
    xr_mss = xr_mss,
    xu_mss = xu_mss
  ))
}
