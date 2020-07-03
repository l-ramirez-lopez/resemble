#' @title Extract predictions from an object of class \code{mbl}
#' @description 
#' 
#' \lifecycle{stable}
#' 
#' Extract predictions from an object of class \code{mbl}
#' @usage
#' get_predictions(object)
#' @param object an object of class \code{mbl} as returned by \code{mbl}
#' @return a data.table of predicted values according to either \code{k} or \code{k_dist}
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @seealso \code{\link{mbl}}
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
get_predictions <- function(object) {
  if (is.na(match("mbl", class(object)))) {
    stop("the object is not of class 'mbl'")
  }

  ext_pred <- function(x, ...) {
    prediction <- x$pred
    return(prediction)
  }
  predictions <- data.table(sapply(object$results, ext_pred))
  return(predictions)
}
