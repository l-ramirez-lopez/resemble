#' @title Evaluate dissimilarity matrices
#'
#' @description
#' \loadmathjax
#' `r lifecycle::badge("stable")`
#'
#' Evaluates a dissimilarity matrix by comparing each observation to its
#' nearest neighbor based on side information. For continuous variables,
#' RMSD and correlation are computed; for categorical variables, the kappa
#' index is used.
#' 
#' @usage
#' diss_evaluate(diss, side_info)
#' 
#' @param diss A symmetric dissimilarity matrix. Alternatively, a vector
#'   containing the lower triangle values (as returned by \code{\link[stats]{dist}}).
#' @param side_info A matrix of side information corresponding to the
#'   observations. Can be numeric (one or more columns) or character (single
#'   column for categorical data).
#'
#' @details
#' This function assesses whether a dissimilarity matrix captures meaningful
#' structure by examining the side information of nearest neighbor pairs
#' (Ramirez-Lopez et al., 2013). If observations that are similar in the
#' dissimilarity space also have similar side information values, the
#' dissimilarity is considered effective.
#'
#' For numeric \code{side_info}, the root mean square of differences (RMSD)
#' between each observation and its nearest neighbor is computed:
#'
#' \mjdeqn{j(i) = NN(x_i, X^{\{-i\}})}{j(i) = NN(x_i, X^{-i})}
#' \mjdeqn{RMSD = \sqrt{\frac{1}{m} \sum_{i=1}^{m} (y_i - y_{j(i)})^2}}{RMSD = sqrt(1/m sum (y_i - y_j(i))^2)}
#'
#' where \mjeqn{NN(x_i, X^{-i})}{NN(x_i, X^{-i})} returns the index of the
#' nearest neighbor of observation \mjeqn{i}{i} (excluding itself),
#' \mjeqn{y_i}{y_i} is the side information value for observation \mjeqn{i}{i},
#' and \mjeqn{m}{m} is the number of observations.
#'
#' For categorical \code{side_info}, the kappa index is computed:
#'
#' \mjdeqn{\kappa = \frac{p_o - p_e}{1 - p_e}}{kappa = (p_o - p_e) / (1 - p_e)}
#'
#' where \mjeqn{p_o}{p_o} is the observed agreement and \mjeqn{p_e}{p_e} is the
#' agreement expected by chance.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{eval}{For numeric side information: a matrix with columns
#'     \code{rmsd} and \code{r} (correlation). For categorical: a matrix
#'     with column \code{kappa}.}
#'   \item{global_eval}{If multiple numeric side information variables are
#'     provided, summary statistics across variables.}
#'   \item{first_nn}{A matrix with the original side information and the
#'     side information of each observation's nearest neighbor.}
#' }
#'
#' @seealso \code{\link{dissimilarity}}, \code{\link{ncomp_by_opc}}
#'
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196,
#' 268-279.
#'
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
#' J.A.M., Scholten, T. 2013b. Distance and similarity-search metrics for use
#' with soil vis-NIR spectra. Geoderma 199, 43-53.
#'
#' @author
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @examples
#' \donttest{
#' library(prospectr)
#' data(NIRsoil)
#'
#' sg <- savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
#' NIRsoil$spc <- sg
#'
#' Yr <- NIRsoil$Nt[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
#'
#' # Compute PCA-based dissimilarity
#' d <- dissimilarity(Xr, diss_method = diss_pca(ncomp = 8))
#'
#' # Evaluate using side information
#' ev <- diss_evaluate(d$dissimilarity, side_info = as.matrix(Yr))
#' ev$eval
#'
#' # Evaluate with multiple side information variables
#' Yr_2 <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' ev_2 <- diss_evaluate(d$dissimilarity, side_info = cbind(Yr, Yr_2))
#' ev_2$eval
#' ev_2$global_eval
#' }
#'
#' @export
diss_evaluate <- function(diss, side_info) {
  # ---------------------------------------------------------------------------
  # Validate side_info
  # ---------------------------------------------------------------------------
  if (!is.matrix(side_info)) {
    stop("'side_info' must be a matrix.", call. = FALSE)
  }
  
  if (is.character(side_info)) {
    if (ncol(side_info) > 1L) {
      stop(
        "For categorical 'side_info', only one column is allowed.",
        call. = FALSE
      )
    }
    side_info <- as.factor(side_info)
    if (nlevels(side_info) < 2L) {
      stop(
        "'side_info' must have at least two categories.",
        call. = FALSE
      )
    }
    get_eval <- .eval_categorical
  } else {
    get_eval <- .eval_continuous
  }
  
  if (any(colSums(!is.na(side_info)) < 4L)) {
    stop(
      "Each 'side_info' variable must have at least 4 non-missing values.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Validate diss
  # ---------------------------------------------------------------------------
  if (!is.vector(diss) && !is.matrix(diss)) {
    stop("'diss' must be a matrix or a vector.", call. = FALSE)
  }
  
  if (is.vector(diss)) {
    expected_len <- (nrow(side_info)^2 - nrow(side_info)) / 2
    if (length(diss) != expected_len) {
      stop(
        "Length of 'diss' does not match a lower triangular matrix ",
        "for ", nrow(side_info), " observations.",
        call. = FALSE
      )
    }
    find_nearest <- which_min_vector
  } else {
    if (nrow(diss) != ncol(diss)) {
      stop("'diss' matrix must be square.", call. = FALSE)
    }
    if (nrow(diss) != nrow(side_info)) {
      stop(
        "Number of rows in 'diss' must match 'side_info'.",
        call. = FALSE
      )
    }
    find_nearest <- which_min
  }
  
  # ---------------------------------------------------------------------------
  # Find nearest neighbors and evaluate
  # ---------------------------------------------------------------------------
  nearest_idx <- find_nearest(diss)
  
  cnms <- colnames(side_info)
  if (is.null(cnms)) {
    cnms <- paste0("side_info_", seq_len(ncol(side_info)))
  }
  eval_result <- get_eval(y = side_info, indices_closest = nearest_idx)
  rownames(eval_result) <- cnms
  
  # ---------------------------------------------------------------------------
  # Build output
  # ---------------------------------------------------------------------------
  first_nn <- side_info[nearest_idx, , drop = FALSE]
  colnames(first_nn) <- paste0("nn_", cnms)
  
  if (ncol(side_info) > 1L) {
    sds <- get_column_sds(side_info[complete.cases(side_info), , drop = FALSE])
    global_eval <- data.frame(
      mean_standardized_rmsd = mean(eval_result[, "rmsd"] / sds),
      mean_r = mean(eval_result[, "r"])
    )
    
    result <- list(
      eval = eval_result,
      global_eval = global_eval,
      first_nn = cbind(side_info, first_nn)
    )
  } else {
    result <- list(
      eval = eval_result,
      first_nn = cbind(side_info, first_nn)
    )
  }
  
  result
}


#' @rdname diss_evaluate
#' @param d Deprecated. Use \code{diss} in \code{diss_evaluate()} instead.
#' @export
sim_eval <- function(d, side_info) {
  .Deprecated("diss_evaluate")
  diss_evaluate(diss = d, side_info = side_info)
}


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.eval_categorical <- function(y, indices_closest) {
  tab <- as.matrix(table(y[!is.na(y)], y[indices_closest][!is.na(y)]))
  total <- sum(tab)
  tab <- tab / total
  p <- rowSums(tab) %*% t(colSums(tab))
  pra <- sum(diag(tab))
  pre <- sum(diag(p))
  kappa <- t((pra - pre) / (1 - pre))
  colnames(kappa) <- "kappa"
  kappa
}

#' @keywords internal
#' @importFrom stats cor
.eval_continuous <- function(y, indices_closest) {
  eval_one <- function(col_idx) {
    obs <- y[, col_idx]
    nn <- y[indices_closest, col_idx]
    rmsd <- sqrt(mean((nn - obs)^2, na.rm = TRUE))
    if (sd(nn, na.rm = TRUE) == 0 || sd(obs, na.rm = TRUE) == 0) {
      r <- NA_real_
    } else {
      r <- cor(nn, obs, use = "complete.obs")
    }
    c(rmsd = rmsd, r = r)
  }
  
  result <- vapply(
    seq_len(ncol(y)),
    eval_one,
    FUN.VALUE = numeric(2)
  )
  t(result)
}
