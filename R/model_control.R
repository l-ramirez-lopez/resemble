#' @title Control parameters for global model fitting
#' @name model_control
#' @description
#' Specifies cross-validation settings for the \code{\link{model}} function.
#'
#' @usage 
#' 
#' model_control(validation_type = c("lgo", "none"), number = 10L, p = 0.75)
#' 
#' @param validation_type a character string specifying the validation method:
#'   \itemize{
#'     \item{\code{"lgo"}: Leave-group-out cross-validation. At each iteration,
#'       a proportion \code{p} of observations is retained for training and the
#'       remainder is used for validation. This is repeated \code{number} times.}
#'     \item{\code{"none"}: No cross-validation is performed.}
#'   }
#' @param number an integer indicating the number of cross-validation 
#'   iterations. Only used when \code{validation_type = "lgo"}. Default is 10.
#' @param p a numeric value between 0 and 1 indicating the proportion of 
#'   observations to retain for training at each cross-validation iteration. 
#'   Only used when \code{validation_type = "lgo"}. Default is 0.75.
#' 
#' @return A list of class \code{"model_control"} containing the specified 
#'   parameters.
#'
#' @seealso \code{\link{model}}
#'
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#'
#' @examples
#' # Default settings (leave-group-out CV with 10 iterations)
#' model_control()
#' 
#' # No cross-validation
#' model_control(validation_type = "none")
#' 
#' # Custom CV settings
#' model_control(validation_type = "lgo", number = 20, p = 0.80)
#'
#' @export
model_control <- function(
    validation_type = c("lgo", "none"),
    number = 10L,
    p = 0.75
) {
  validation_type <- match.arg(validation_type)
  
  
  if (!is.numeric(number) || length(number) != 1L || number < 1L) {
    stop("'number' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1) {
    stop("'p' must be a numeric value between 0 and 1 (exclusive).", 
         call. = FALSE)
  }
  
  structure(
    list(
      validation_type = validation_type,
      number          = as.integer(number),
      p               = p
    ),
    class = "model_control"
  )
}


#' @noRd
#' @export
print.model_control <- function(x, ...) {
  cat("Model control parameters:\n")
  cat("  validation_type :", x$validation_type, "\n")
  if (x$validation_type != "none") {
    cat("  number          :", x$number, "\n")
    cat("  p               :", x$p, "\n")
  }
  invisible(x)
}


