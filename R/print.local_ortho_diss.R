#' @title Print method for an object of class \code{ortho_diss}
#' @description Prints the content of an object of class \code{ortho_diss}
#' @usage \method{print}{local_ortho_diss}(x, ...)
#' @param x an object of class \code{localortho_diss} (returned by 
#' \code{ortho_diss} when it uses \code{.local = TRUE}).
#' @param ... arguments to be passed to methods (not yet functional).
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @keywords internal
#' @export

######################################################################
# resemble
# Copyright (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
######################################################################

## History:
## 09.03.2014 Leo     The tryMod function was removed

print.local_ortho_diss <- function(x, ...) {
  obj <- x
  if (is.list(obj)) {
    obj$dissimilarity <- as.matrix(round(obj$dissimilarity, getOption("digits")))
    obj$dissimilarity[is.na(obj$dissimilarity)] <- "*"
    dm <- format(obj$dissimilarity, digits = getOption("digits"), justify = "right")
    print(list(n.components = object$n_components, 
               loc_n_components = object$neighborhood_info$local_n_components, 
               dissimilarity = noquote(dm)))
    cat("*: Not a (local) neighbor")
  }
  if (is.matrix(obj)) {
    object <- as.matrix(round(obj, getOption("digits")))
    object[is.na(obj)] <- "*"
    dm <- format(object, digits = getOption("digits"), justify = "right")
    print(dm, quote = FALSE)
    cat("*: not a neighbor")
  }
}

#' @title Print method for an object of class \code{local_ortho_diss}
#' @param x \code{local_ortho_diss} matrix
#' @param rows the indices of the rows
#' @param columns the indices of the columns
#' @param drop drop argument
#' @param ... not used
#' @description prints the subsets of localortho_diss objects
#' @keywords internal
#' @export
"[.local_ortho_diss" <- function(x, rows, columns, drop = FALSE, ...) {
  object <- x
  if (!is.logical(drop)) {
    drop <- FALSE
  }
  class(object) <- NULL
  obj <- object[rows, columns, drop = drop]
  if (!drop) {
    class(obj) <- "localortho_diss"
  }
  return(obj)
}
