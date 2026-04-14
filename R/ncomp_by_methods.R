# =============================================================================
# ncomp_selection constructors
# =============================================================================
# Public:   ncomp_by_var()    — select by minimum per-component variance
#           ncomp_by_cumvar() — select by cumulative variance
#           ncomp_by_opc()    — select by OPC method
#           ncomp_fixed()     — fix number of components manually


#' @title Select number of components by minimum per-component variance
#'
#' @description
#' Retains components that individually explain at least \code{min_var}
#' proportion of variance. Use with \code{diss_pca()} or \code{diss_pls()}.
#'
#' @param min_var A number in (0, 1]. Minimum proportion of variance a single
#'   component must explain to be retained. Default is \code{0.01} (1\%).
#'
#' @return An object of class \code{c("ncomp_by_var", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_cumvar}}, \code{\link{ncomp_by_opc}},
#'   \code{\link{ncomp_fixed}}
#' @examples
#' diss_pca(ncomp = ncomp_by_var(0.01))
#' diss_pls(ncomp = ncomp_by_var(0.05))
#' @export
ncomp_by_var <- function(min_var = 0.01) {
  if (!is.numeric(min_var) || length(min_var) != 1L ||
      min_var <= 0 || min_var > 1) {
    stop("'min_var' must be a single number in (0, 1].")
  }
  structure(
    list(method = "var", value = min_var),
    class = c("ncomp_by_var", "ncomp_selection")
  )
}


#' @title Select number of components by cumulative variance explained
#'
#' @description
#' Retains the minimum number of components whose combined explained variance
#' reaches at least \code{min_cumvar}. Use with \code{diss_pca()} or
#' \code{diss_pls()}.
#'
#' @param min_cumvar A number in (0, 1]. Minimum cumulative proportion of
#'   variance the retained components must collectively explain. Default is
#'   \code{0.99} (99\%).
#'
#' @return An object of class \code{c("ncomp_by_cumvar", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_opc}},
#'   \code{\link{ncomp_fixed}}
#' @examples
#' diss_pca(ncomp = ncomp_by_cumvar(0.99))
#' diss_pls(ncomp = ncomp_by_cumvar(0.95))
#' @export
ncomp_by_cumvar <- function(min_cumvar = 0.99) {
  if (!is.numeric(min_cumvar) || length(min_cumvar) != 1L ||
      min_cumvar <= 0 || min_cumvar > 1) {
    stop("'min_cumvar' must be a single number in (0, 1].")
  }
  structure(
    list(method = "cumvar", value = min_cumvar),
    class = c("ncomp_by_cumvar", "ncomp_selection")
  )
}


#' @title Select number of components using the OPC method
#'
#' @description
#' Selects the optimal number of components based on the side information
#' concept (Ramirez-Lopez et al., 2013). The optimal number is the one for
#' which the dissimilarity matrix minimizes differences between each
#' observation's \code{Yr} value and that of its closest neighbor.
#'
#' Requires \code{Yr} to be passed to \code{dissimilarity()}.
#'
#' @param max_ncomp A positive integer. Maximum number of components to
#'   evaluate. Default is \code{40}.
#'
#' @return An object of class \code{c("ncomp_by_opc", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_cumvar}},
#'   \code{\link{ncomp_fixed}}
#' @examples
#' diss_pca(ncomp = ncomp_by_opc(40))
#' diss_pls(ncomp = ncomp_by_opc(30))
#' @export
ncomp_by_opc <- function(max_ncomp = 40) {
  if (!is.numeric(max_ncomp) || length(max_ncomp) != 1L) {
    stop("'max_ncomp' must be a single positive integer.")
  }
  max_ncomp <- as.integer(max_ncomp)
  if (max_ncomp < 1L) {
    stop("'max_ncomp' must be at least 1.")
  }
  structure(
    list(method = "opc", value = max_ncomp),
    class = c("ncomp_by_opc", "ncomp_selection")
  )
}


#' @title Fix the number of components manually
#'
#' @description
#' Retains a fixed, user-specified number of components. Use when you know
#' exactly how many components are appropriate for your data.
#'
#' @param ncomp A positive integer. Number of components to retain.
#'
#' @return An object of class \code{c("ncomp_fixed", "ncomp_selection")}.
#' @seealso \code{\link{ncomp_by_var}}, \code{\link{ncomp_by_cumvar}},
#'   \code{\link{ncomp_by_opc}}
#' @examples
#' diss_pca(ncomp = ncomp_fixed(10))
#' diss_pls(ncomp = ncomp_fixed(5))
#' @export
ncomp_fixed <- function(ncomp) {
  if (missing(ncomp)) {
    stop("'ncomp' is required for ncomp_fixed().")
  }
  if (!is.numeric(ncomp) || length(ncomp) != 1L) {
    stop("'ncomp' must be a single positive integer.")
  }
  ncomp <- as.integer(ncomp)
  if (ncomp < 1L) {
    stop("'ncomp' must be at least 1.")
  }
  structure(
    list(method = "manual", value = ncomp),
    class = c("ncomp_fixed", "ncomp_selection")
  )
}


#' @export
print.ncomp_selection <- function(x, ...) {
  label <- switch(
    class(x)[[1]],
    ncomp_by_var    = paste("by per-component variance >=", x$value),
    ncomp_by_cumvar = paste("by cumulative variance >=",    x$value),
    ncomp_by_opc    = paste("by OPC (max components:",      x$value, ")"),
    ncomp_fixed     = paste("fixed at",                     x$value, "component(s)")
  )
  cat("Component selection:", label, "\n")
  invisible(x)
}

