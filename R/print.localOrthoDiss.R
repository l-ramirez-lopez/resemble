#' @title Print method for an object of class \code{orthoDiss}
#' @aliases print.localOrthoDiss
#' @usage \method{print}{localOrthoDiss}(x, ...)
#' @param x an object of class \code{localOrthoDiss} (returned by \code{orthoDiss} when it uses \code{local = TRUE}). 
#' @param ... arguments to be passed to methods (not yet functional).
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @export

print.localOrthoDiss <- function(x,...){
  obj <- x
  if(is.list(obj))
  {
    obj$dissimilarity <- as.matrix(round(obj$dissimilarity, getOption("digits")))
    obj$dissimilarity[is.na(obj$dissimilarity)] <- "*"
    dm <- format(obj$dissimilarity, digits = getOption("digits"), justify = "right")
    print(list(n.components = object$n.components, loc.n.components = obj$loc.n.components, dissimilarity = noquote(dm)))
    cat("*: local non-neighbor sample")
  }
  if(is.matrix(obj))
  {
    object <- as.matrix(round(obj, getOption("digits")))
    object[is.na(obj)] <- "*"
    dm <- format(object, digits = getOption("digits"), justify = "right")
    print(dm, quote = FALSE)
    cat("*: local non-neighbor sample")
  }
}

"[.localOrthoDiss" <- function(x, rr, cl, drop = FALSE, ...){
  object <- x
  if(!is.logical(drop))
    drop <- FALSE
  class(object) <- NULL
  obj <- tryMod(object[rr,cl, drop = drop])
  if(!drop)
    class(obj) <- "localOrthoDiss"
  return(obj)
}  

tryMod <- function (expr, addMss = NULL, silent = FALSE) 
{
  tryCatch(expr, error = function(e) {
    call <- conditionCall(e)
    prefix <- "Error: "
    msg <- paste0(prefix, conditionMessage(e), "\n")
    .Internal(seterrmessage(msg[1L]))
    if (!silent && identical(getOption("show.error.messages"), TRUE)) 
    {
      cat(msg, file = stderr())
      .Internal(printDeferredWarnings())
    }
    invisible(structure(msg, class = "try-error", condition = e))
  })
}