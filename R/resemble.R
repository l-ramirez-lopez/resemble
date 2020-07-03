#' This is the version \code{1.3} (\code{'milagros'}) of the package. It 
#' implements a number of \code{R} functions useful for 
#' modeling complex spectral spectra (e.g. NIR, IR). 
#' The package includes functions for projecting spectral data 
#' onto orthogonal spaces, computing spectral dissimilarity
#' matrices, removing irrelevant spectra from a reference set, 
#' and modeling spectral data using memory-based learning.
#' 
#' The functions available for projecting the spectra are:
#' \itemize{
#'   \item{\code{\link{orthoProjection}}} 
#'   \item{\code{\link{pcProjection}}} 
#'   \item{\code{\link{plsProjection}}} 
#'   \item{\code{\link{predict.orthoProjection}}} 
#'   }
#' The functions available for computing similarity/dissimilarity 
#' matrices are:
#' \itemize{
#'   \item{\code{\link{fDiss}}} 
#'   \item{\code{\link{corDiss}}} 
#'   \item{\code{\link{sid}}} 
#'   \item{\code{\link{orthoDiss}}} 
#'   }
#' The functions available for evaluating similarity/dissimilarity 
#' matrices are:
#' \itemize{
#'   \item{\code{\link{simEval}}} 
#'   }
#' The functions available for removing irrelevant spectra from a 
#' reference set are:
#' \itemize{
#'   \item{\code{\link{neigCleaning}}} 
#'   }
#' The functions available for modeling spectral data:
#' \itemize{
#'   \item{\code{\link{mblControl}}}
#'   \item{\code{\link{mbl}}} 
#'   \item{\code{\link{rs_control}}} 
#'   \item{\code{\link{rslocal}}} 
#'   \item{\code{\link{predict.rslocal}}} 
#'   }
#' Other supplementary functions are:
#' \itemize{
#'   \item{\code{\link{plot.mbl}}}
#'   \item{\code{\link{plot.orthoProjection}}}
#'   }
#' @docType package
#' @name resemble-package
#' @aliases resemble-package resemble
#' @title Overview of the functions in the resemble package
#' @import Rcpp RcppArmadillo foreach iterators
#' @useDynLib resemble
#' @author 
#' Leonardo Ramirez-Lopez \email{ramirez.lopez.leo@@gmail.com}, 
#' Antoine Stevens, 
#' Craig Lobsey, 
#' Raphael Viscarra-Rossel

######################################################################
# resemble
# Copyrigth (C) 2020 Leonardo Ramirez-Lopez
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


NULL
