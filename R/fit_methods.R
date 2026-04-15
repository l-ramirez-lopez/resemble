#' @title Local fitting method constructors
#' @name fit_methods
#' @aliases fit_pls
#' @aliases fit_wapls
#' @aliases fit_gpr
#' @description
#' 
#' \loadmathjax
#' 
#' These functions create configuration objects that specify how local 
#' regression models are fitted within the \code{\link{mbl}} function.
#' @usage
#' fit_pls(ncomp, method = c("pls", "mpls", "simpls"), 
#'         scale = FALSE, max_iter = 100L, tol = 1e-6)
#'
#' fit_wapls(min_ncomp, max_ncomp, method = c("mpls", "pls", "simpls"),
#'           scale = FALSE, max_iter = 100L, tol = 1e-6)
#'
#' fit_gpr(noise_variance = 0.001, center = TRUE, scale = TRUE)
#'
#' @param ncomp an integer indicating the number of PLS components to use
#' in local regressions when \code{fit_pls} is used.
#' @param min_ncomp an integer indicating the minimum number of PLS components
#' to use in local regressions when \code{fit_wapls} is used. See details.
#' @param max_ncomp an integer indicating the maximum number of PLS components
#' to use in local regressions when \code{fit_wapls} is used. See details.
#' @param method a character string indicating the PLS algorithm to use. 
#' Options are:
#' \itemize{
#'   \item{\code{'pls'}: standard PLS using covariance between X and Y for 
#'     weight computation (NIPALS algorithm).}
#'   \item{\code{'mpls'}: modified PLS using correlation between X and Y for 
#'     weight computation (NIPALS algorithm). See Shenk and Westerhaus (1991).}
#'   \item{\code{'simpls'}: SIMPLS algorithm (de Jong, 1993). Computationally 
#'     faster as it avoids iterative X deflation. Parameters \code{max_iter} 
#'     and \code{tol} are ignored when this method is used.}
#' }
#' Default is \code{'pls'} for \code{fit_pls} and \code{'mpls'} for 
#' \code{fit_wapls}.
#' @param scale logical indicating whether predictors must be scaled. 
#' Default is \code{FALSE} for PLS methods and \code{TRUE} for GPR.
#' @param max_iter an integer indicating the maximum number of iterations 
#' for convergence in the NIPALS algorithm. Only used when 
#' \code{method = 'pls'} or \code{method = 'mpls'}. Default is 100.
#' @param tol a numeric value indicating the convergence tolerance for 
#' calculating scores in the NIPALS algorithm. Only used when 
#' \code{method = 'pls'} or \code{method = 'mpls'}. Default is 1e-6.
#' @param noise_variance a numeric value indicating the variance of the noise
#' for Gaussian process local regressions (\code{fit_gpr}). Default is 0.001.
#' @param center logical indicating whether predictors should be centered 
#' before fitting. Only used for \code{fit_gpr}. Default is \code{TRUE}.
#'
#' @details
#' These functions create configuration objects that are passed to 
#' \code{\link{mbl}} to specify how local regression models are fitted.
#'
#' There are three fitting methods available:
#' 
#' \subsection{Partial least squares (\code{fit_pls})}{
#' Uses orthogonal scores partial least squares regression. Three algorithm 
#' variants are available:
#' \itemize{
#'   \item{\strong{Standard PLS} (\code{method = 'pls'}): Uses the NIPALS 
#'     algorithm with covariance-based weights.}
#'   \item{\strong{Modified PLS} (\code{method = 'mpls'}): Uses the NIPALS 
#'     algorithm with correlation-based weights. Proposed by Shenk and 
#'     Westerhaus (1991), this approach gives equal influence to all 
#'     predictors regardless of their variance scale.}
#'   \item{\strong{SIMPLS} (\code{method = 'simpls'}): Uses the SIMPLS 
#'     algorithm (de Jong, 1993), which deflates the cross-product matrix 
#'     rather than X itself. This is computationally faster, especially for 
#'     wide matrices, and produces identical predictions to standard PLS.}
#' }
#' The only parameter to optimise is the number of PLS components 
#' (\code{ncomp}).
#' }
#' 
#' \subsection{Weighted average PLS (\code{fit_wapls})}{
#' This method was developed by Shenk et al. (1997) and is used as the 
#' regression method in the LOCAL algorithm. It fits multiple PLS models 
#' using different numbers of components (from \code{min_ncomp} to 
#' \code{max_ncomp}). The final prediction is a weighted average of 
#' predictions from all models, where the weight for component \mjeqn{j}{j} 
#' is:
#'
#' \mjdeqn{w_{j} = \frac{1}{s_{1:j} \times g_{j}}}{w_j = 1/(s_{1:j} * g_j)}
#'
#' where \mjeqn{s_{1:j}}{s_{1:j}} is the root mean square of the spectral 
#' reconstruction error of the target observation(s) when \mjeqn{j}{j} PLS 
#' components are used, and \mjeqn{g_{j}}{g_j} is the root mean square of 
#' the squared regression coefficients for the \mjeqn{j}{j}th component.
#'
#' The same algorithm variants (\code{'pls'}, \code{'mpls'}, \code{'simpls'}) 
#' are available. The default is \code{'mpls'} following the original LOCAL 
#' implementation.
#' }
#' 
#' \subsection{Gaussian process regression (\code{fit_gpr})}{
#' Gaussian process regression is a non-parametric Bayesian method 
#' characterised by a mean and covariance function. This implementation uses 
#' a dot product covariance.
#'
#' The prediction vector \mjeqn{A}{A} is computed from training data 
#' (\mjeqn{X}{X}, \mjeqn{Y}{Y}) as:
#'
#' \mjdeqn{A = (X X^{T} + \sigma^2 I)^{-1} Y}{A = (X X^T + sigma^2 I)^{-1} Y}
#'
#' where \mjeqn{\sigma^2}{sigma^2} is the noise variance and \mjeqn{I}{I} is 
#' the identity matrix. Prediction for a new observation 
#' \mjeqn{x_{u}}{x_u} is:
#'
#' \mjdeqn{\hat{y}_{u} = x_{u} X^{T} A}{hat y_u = x_u X^T A}
#'
#' The only parameter is the noise variance (\code{noise_variance}).
#' }
#'
#' @return An object of class \code{c("fit_<method>", "fit_method")} 
#' containing the specified parameters. This object is passed to 
#' \code{\link{mbl}} to configure local model fitting.
#' 
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' 
#' @references
#' de Jong, S. (1993). SIMPLS: An alternative approach to partial least 
#' squares regression. Chemometrics and Intelligent Laboratory Systems, 
#' 18(3), 251-263.
#' 
#' Rasmussen, C.E., Williams, C.K. (2006). Gaussian Processes for Machine 
#' Learning. MIT Press.
#' 
#' Shenk, J.S., & Westerhaus, M.O. (1991). Populations structuring of near 
#' infrared spectra and modified partial least squares regression. Crop 
#' Science, 31(6), 1548-1555.
#'
#' Shenk, J., Westerhaus, M., & Berzaghi, P. (1997). Investigation of a LOCAL
#' calibration procedure for near infrared instruments. Journal of Near 
#' Infrared Spectroscopy, 5, 223-232.
#' 
#' Westerhaus, M. (2014). Eastern Analytical Symposium Award for outstanding 
#' achievements in near infrared spectroscopy: my contributions to near 
#' infrared spectroscopy. NIR news, 25(8), 16-20.
#' 
#' @seealso \code{\link{mbl}}
#' 
#' @examples
#' # PLS with 10 components using standard algorithm
#' fit_pls(ncomp = 10)
#' 
#' # PLS with modified algorithm (correlation-based weights)
#' fit_pls(ncomp = 10, method = "mpls")
#' 
#' # PLS with SIMPLS (faster, no iteration)
#' fit_pls(ncomp = 10, method = "simpls")
#' 
#' # Weighted average PLS (LOCAL-style)
#' fit_wapls(min_ncomp = 3, max_ncomp = 12)
#' 
#' # Weighted average PLS with SIMPLS
#' fit_wapls(min_ncomp = 3, max_ncomp = 15, method = "simpls")
#' 
#' # Gaussian process regression
#' fit_gpr()
#' fit_gpr(noise_variance = 0.01)
#' 
NULL


# =============================================================================
# Fitting method constructors
# =============================================================================
# Public:
#   fit_pls()
#   fit_wapls()
#   fit_gpr()
#
# Internal:
#   .new_fit_method()    - shared constructor
#   print.fit_method     - shared print (fallback)
#   print.fit_pls
#   print.fit_wapls
#   print.fit_gpr


# -----------------------------------------------------------------------------
# Internal constructor
# -----------------------------------------------------------------------------
.new_fit_method <- function(fit_method, ...) {
  
  structure(
    list(fit_method = fit_method, ...),
    class = c(paste0("fit_", fit_method), "fit_method")
  )
}


# -----------------------------------------------------------------------------
# fit_pls
# -----------------------------------------------------------------------------
#' @export
fit_pls <- function(
    ncomp,
    method = c("pls", "mpls", "simpls"),
    scale = FALSE,
    max_iter = 100L,
    tol = 1e-6
) {
  if (missing(ncomp)) {
    stop("'ncomp' is required.", call. = FALSE)
  }
  
  if (!is.numeric(ncomp) || length(ncomp) != 1L || ncomp < 1L) {
    stop("'ncomp' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.logical(scale) || length(scale) != 1L) {
    stop("'scale' must be TRUE or FALSE.", call. = FALSE)
  }
  
  method <- match.arg(method)
  
  .new_fit_method(
    "pls",
    ncomp    = as.integer(ncomp),
    method   = method,
    scale    = scale,
    max_iter = as.integer(max_iter),
    tol      = tol
  )
}

#' @noRd
#' @export
print.fit_pls <- function(x, ...) {
  cat(
    "Fitting method: pls\n",
    "  ncomp    :", x$ncomp, "\n",
    "  method   :", x$method, "\n",
    "  scale    :", x$scale, "\n",
    "  max_iter :", x$max_iter, "\n",
    "  tol      :", x$tol, "\n"
  )
  invisible(x)
}


# -----------------------------------------------------------------------------
# fit_wapls
# -----------------------------------------------------------------------------
#' @export
fit_wapls <- function(
    min_ncomp,
    max_ncomp,
    method = c("mpls", "pls", "simpls"),
    scale = FALSE,
    max_iter = 100L,
    tol = 1e-6
) {
  if (missing(min_ncomp) || missing(max_ncomp)) {
    stop("Both 'min_ncomp' and 'max_ncomp' are required.", call. = FALSE)
  }
  
  if (!is.numeric(min_ncomp) || length(min_ncomp) != 1L || min_ncomp < 1L) {
    stop("'min_ncomp' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.numeric(max_ncomp) || length(max_ncomp) != 1L || max_ncomp < 1L) {
    stop("'max_ncomp' must be a positive integer.", call. = FALSE)
  }
  
  if (min_ncomp >= max_ncomp) {
    stop("'min_ncomp' must be less than 'max_ncomp'.", call. = FALSE)
  }
  
  if (!is.logical(scale) || length(scale) != 1L) {
    stop("'scale' must be TRUE or FALSE.", call. = FALSE)
  }
  
  method <- match.arg(method)
  
  .new_fit_method(
    "wapls",
    min_ncomp = as.integer(min_ncomp),
    max_ncomp = as.integer(max_ncomp),
    method    = method,
    scale     = scale,
    max_iter  = as.integer(max_iter),
    tol       = tol
  )
}

#' @noRd
#' @export
print.fit_wapls <- function(x, ...) {
  cat(
    "Fitting method: wapls\n",
    "  min_ncomp :", x$min_ncomp, "\n",
    "  max_ncomp :", x$max_ncomp, "\n",
    "  method    :", x$method, "\n",
    "  scale     :", x$scale, "\n",
    "  max_iter  :", x$max_iter, "\n",
    "  tol       :", x$tol, "\n"
  )
  invisible(x)
}



# -----------------------------------------------------------------------------
# fit_gpr
# -----------------------------------------------------------------------------

#' @export
fit_gpr <- function(
    noise_variance = 0.001,
    center = TRUE,
    scale = TRUE
) {
  if (!is.numeric(noise_variance) || length(noise_variance) != 1L) {
    stop("'noise_variance' must be a single numeric value.", call. = FALSE)
  }
  
  if (noise_variance <= 0) {
    stop("'noise_variance' must be positive.", call. = FALSE)
  }
  
  if (!is.logical(center) || length(center) != 1L) {
    stop("'center' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(scale) || length(scale) != 1L) {
    stop("'scale' must be TRUE or FALSE.", call. = FALSE)
  }
  
  .new_fit_method(
    "gpr",
    noise_variance = noise_variance,
    center         = center,
    scale          = scale
  )
}

#' @noRd
#' @export
print.fit_gpr <- function(x, ...) {
  cat(
    "Fitting method: gpr\n",
    "  noise_variance :", x$noise_variance, "\n",
    "  center         :", x$center, "\n",
    "  scale          :", x$scale, "\n"
  )
  invisible(x)
}


# -----------------------------------------------------------------------------
# Fallback print for fit_method
# -----------------------------------------------------------------------------
#' @noRd
#' @export
print.fit_method <- function(x, ...) {
  cat("Fitting method:", x$fit_method, "\n")
  for (nm in setdiff(names(x), "fit_method")) {
    cat(" ", nm, ":", x[[nm]], "\n")
  }
  invisible(x)
}
