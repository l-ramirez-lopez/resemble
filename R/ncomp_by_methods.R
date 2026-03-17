# =============================================================================
# ncomp_selection constructors — revised
# =============================================================================

#' @title Select number of components by minimum per-component variance
#'
#' @description
#' Retains components that individually explain at least \code{min_var}
#' proportion of variance, up to a maximum of \code{max_ncomp} components.
#'
#' @param min_var Numeric in (0, 1]. Minimum variance a component must
#'   explain to be retained. Default \code{0.01}.
#' @param max_ncomp Positive integer. Maximum number of components to
#'   compute/evaluate. Default \code{40}.
#'
#' @return An object of class \code{c("ncomp_by_var", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_cumvar}}, \code{\link{ncomp_by_opc}},
#'   \code{\link{ncomp_fixed}}
#' @examples
#' diss_pca(ncomp = ncomp_by_var(0.01))
#' diss_pca(ncomp = ncomp_by_var(0.05, max_ncomp = 20))
#' @export
ncomp_by_var <- function(min_var = 0.01, max_ncomp = 40L) {
  if (!is.numeric(min_var) || length(min_var) != 1L ||
      min_var <= 0 || min_var > 1) {
    stop("'min_var' must be a single number in (0, 1].")
  }
  .validate_max_ncomp(max_ncomp)
  
  structure(
    list(
      method = "var",
      min_var = min_var,
      max_ncomp = as.integer(max_ncomp)
    ),
    class = c("ncomp_by_var", "ncomp_selection")
  )
}


#' @title Select number of components by cumulative variance explained
#'
#' @description
#' Retains the minimum number of components whose combined explained variance
#' reaches at least \code{min_cumvar}, up to \code{max_ncomp} components.
#'
#' @param min_cumvar Numeric in (0, 1]. Minimum cumulative variance.
#'   Default \code{0.99}.
#' @param max_ncomp Positive integer. Maximum components to compute.
#'   Default \code{40}.
#'
#' @return An object of class \code{c("ncomp_by_cumvar", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_opc}},
#'   \code{\link{ncomp_fixed}}
#' @examples
#' diss_pca(ncomp = ncomp_by_cumvar(0.99))
#' diss_pca(ncomp = ncomp_by_cumvar(0.95, max_ncomp = 30))
#' @export
ncomp_by_cumvar <- function(min_cumvar = 0.99, max_ncomp = 40L) {
  if (!is.numeric(min_cumvar) || length(min_cumvar) != 1L ||
      min_cumvar <= 0 || min_cumvar > 1) {
    stop("'min_cumvar' must be a single number in (0, 1].")
  }
  .validate_max_ncomp(max_ncomp)
  
  structure(
    list(
      method = "cumvar",
      min_cumvar = min_cumvar,
      max_ncomp = as.integer(max_ncomp)
    ),
    class = c("ncomp_by_cumvar", "ncomp_selection")
  )
}


#' @title Select number of components using the OPC method
#'
#' @description
#' Optimized principal component selection based on side information
#' (Ramirez-Lopez et al., 2013). Requires \code{Yr} in \code{dissimilarity()}.
#'
#' @param max_ncomp Positive integer. Maximum components to evaluate.
#'   Default \code{40}.
#'
#' @return An object of class \code{c("ncomp_by_opc", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_cumvar}},
#'   \code{\link{ncomp_fixed}}
#' @examples
#' diss_pca(ncomp = ncomp_by_opc(40))
#' diss_pls(ncomp = ncomp_by_opc(30))
#' @export
ncomp_by_opc <- function(max_ncomp = 40L) {
  .validate_max_ncomp(max_ncomp)
  
  structure(
    list(
      method = "opc",
      max_ncomp = as.integer(max_ncomp)
    ),
    class = c("ncomp_by_opc", "ncomp_selection")
  )
}


#' @title Fix the number of components
#'
#' @description
#' Use exactly \code{ncomp} components. No automatic selection.
#'
#' @param ncomp Positive integer. Exact number of components.
#'
#' @return An object of class \code{c("ncomp_fixed", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_cumvar}},
#'   \code{\link{ncomp_by_opc}}
#' @examples
#' diss_pca(ncomp = ncomp_fixed(10))
#' diss_pca(ncomp = 10)
#' @export
ncomp_fixed <- function(ncomp) {
  if (missing(ncomp)) {
    stop("'ncomp' is required.")
  }
  if (!is.numeric(ncomp) || length(ncomp) != 1L || ncomp < 1L) {
    stop("'ncomp' must be a single positive integer.")
  }
  
  structure(
    list(
      method = "manual",
      ncomp = as.integer(ncomp)
    ),
    class = c("ncomp_fixed", "ncomp_selection")
  )
}


# --- internal validator ------------------------------------------------------
.validate_max_ncomp <- function(max_ncomp) {
  if (!is.numeric(max_ncomp) || length(max_ncomp) != 1L || max_ncomp < 1L) {
    stop("'max_ncomp' must be a single positive integer.")
  }
  invisible(NULL)
}


# --- print method ------------------------------------------------------------
#' @export
print.ncomp_selection <- function(x, ...) {
  label <- switch(
    class(x)[[1]],
    ncomp_by_var = sprintf(
      "by per-component variance >= %.2f (max: %d)",
      x$min_var, x$max_ncomp
    ),
    ncomp_by_cumvar = sprintf(
      "by cumulative variance >= %.2f (max: %d)",
      x$min_cumvar, x$max_ncomp
    ),
    ncomp_by_opc = sprintf("by OPC (max: %d)", x$max_ncomp),
    ncomp_fixed = sprintf("fixed: %d", x$ncomp)
  )
  cat("Component selection:", label, "\n")
  invisible(x)
}