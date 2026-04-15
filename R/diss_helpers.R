# -----------------------------------------------------------------------------
# Shared constructor helper — not exported
# Avoids repeating the same center/scale/precision validation three times
# -----------------------------------------------------------------------------
.new_diss_method <- function(method_name, center, scale, extra = list()) {
  
  if (!is.logical(center) || length(center) != 1L) {
    stop("'center' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(scale) || length(scale) != 1L) {
    stop("'scale' must be a single logical value (TRUE or FALSE).")
  }
  
  structure(
    c(list(method = method_name, center = center, scale = scale), extra),
    class = c(paste0("diss_", method_name), "diss_method")
  )
}

# -----------------------------------------------------------------------------
# Shared print method for all three — dispatches on "diss_method" base class
# -----------------------------------------------------------------------------
#' @export
print.diss_method <- function(x, ...) {
  cat(
    "Dissimilarity method:", x$method, "\n",
    " center :", x$center, "\n",
    " scale  :", x$scale, "\n"
  )
  invisible(x)
}
