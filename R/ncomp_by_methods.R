#' @title Component selection methods
#'
#' @description
#' Constructor functions for specifying how to select the number of components
#' in projection-based dissimilarity methods ([diss_pca()], [diss_pls()]).
#'
#' @param min_var Numeric in (0, 1]. Minimum variance a single component must
#'   explain to be retained.
#' @param min_cumvar Numeric in (0, 1]. Minimum cumulative variance that the
#'   retained components must explain.
#' @param max_ncomp Positive integer. Maximum number of components to
#'   compute or evaluate.
#' @param ncomp Positive integer. Exact number of components to use.
#'
#' @details
#' Four selection methods are available:
#'
#' \describe{
#'   \item{`ncomp_by_var()`}{Retains components that individually explain at
#'     least `min_var` proportion of variance.}
#'   \item{`ncomp_by_cumvar()`}{Retains the minimum number of components whose
#'     combined explained variance reaches `min_cumvar`.
#'   }
#'   \item{`ncomp_by_opc()`}{Optimized principal component selection based on
#'     side information (Ramirez-Lopez et al., 2013). The optimal number of
#'     components minimizes the RMSD between each observation's response and
#'     its nearest neighbor's response in the projected space. Requires `Yr`.}
#'   \item{`ncomp_fixed()`}{Uses exactly `ncomp` components with no automatic
#'     selection. Equivalent to passing an integer directly.}
#' }
#'
#' At runtime, `max_ncomp` is capped at `min(max_ncomp, nrow(X), ncol(X))`.
#'
#' @return
#' An object of class `"ncomp_selection"` with a subclass indicating the
#' method:
#' \itemize{
#'   \item `ncomp_by_var`: class `c("ncomp_by_var", "ncomp_selection")`
#'   \item `ncomp_by_cumvar`: class `c("ncomp_by_cumvar", "ncomp_selection")`
#'   \item `ncomp_by_opc`: class `c("ncomp_by_opc", "ncomp_selection")`
#'   \item `ncomp_fixed`: class `c("ncomp_fixed", "ncomp_selection")`
#' }
#'
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196,
#' 268-279.
#'
#' @seealso [diss_pca()], [diss_pls()], [dissimilarity()]
#'
#' @examples
#' # Retain components explaining >= 1% variance each
#' ncomp_by_var(0.01)
#'
#' # Retain enough components for 99% cumulative variance
#' ncomp_by_cumvar(0.99)
#'
#' # Optimize using side information (requires Yr)
#' ncomp_by_opc(max_ncomp = 40)
#'
#' # Fix at exactly 10 components
#' ncomp_fixed(10)
#'
#' # Usage in dissimilarity constructors
#' diss_pca(ncomp = ncomp_by_var(0.01))
#' diss_pca(ncomp = ncomp_by_opc())
#' diss_pca(ncomp = 10)  
#'
#' @name ncomp_selection
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
NULL

#' @rdname ncomp_selection
#' @export
ncomp_by_var <- function(min_var = 0.01, max_ncomp = 40L) {
  
  if (!is.numeric(min_var) || length(min_var) != 1L ||
      min_var <= 0 || min_var > 1) {
    stop("'min_var' must be a single number in (0, 1].", call. = FALSE)
  }
  .validate_max_ncomp(max_ncomp)
  
  structure(
    list(
      min_var = min_var,
      max_ncomp = as.integer(max_ncomp)
    ),
    class = c("ncomp_by_var", "ncomp_selection")
  )
}

#' @rdname ncomp_selection
#' @export
ncomp_by_cumvar <- function(min_cumvar = 0.99, max_ncomp = 40L) {
  if (!is.numeric(min_cumvar) || length(min_cumvar) != 1L ||
      min_cumvar <= 0 || min_cumvar > 1) {
    stop("'min_cumvar' must be a single number in (0, 1].", call. = FALSE)
  }
  .validate_max_ncomp(max_ncomp)
  
  structure(
    list(
      min_cumvar = min_cumvar,
      max_ncomp = as.integer(max_ncomp)
    ),
    class = c("ncomp_by_cumvar", "ncomp_selection")
  )
}

#' @rdname ncomp_selection
#' @export
ncomp_by_opc <- function(max_ncomp = 40L) {
  .validate_max_ncomp(max_ncomp)
  
  structure(
    list(max_ncomp = as.integer(max_ncomp)),
    class = c("ncomp_by_opc", "ncomp_selection")
  )
}

#' @rdname ncomp_selection
#' @export
ncomp_fixed <- function(ncomp) {
  if (missing(ncomp)) {
    stop("'ncomp' is required.", call. = FALSE)
  }
  if (!is.numeric(ncomp) || length(ncomp) != 1L || ncomp < 1L) {
    stop("'ncomp' must be a single positive integer.", call. = FALSE)
  }
  
  structure(
    list(ncomp = as.integer(ncomp)),
    class = c("ncomp_fixed", "ncomp_selection")
  )
}

# --- internal validator ------------------------------------------------------
.validate_max_ncomp <- function(max_ncomp) {
  if (!is.numeric(max_ncomp) || length(max_ncomp) != 1L || max_ncomp < 1L) {
    stop("'max_ncomp' must be a single positive integer.", call. = FALSE)
  }
  invisible(NULL)
}

# --- print method ------------------------------------------------------------
#' @export
print.ncomp_selection <- function(x, ...) {
  label <- switch(
    class(x)[[1L]],
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
