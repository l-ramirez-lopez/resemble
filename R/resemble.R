#' @useDynLib resemble
#' @import lifecycle
#' @import Rcpp
#' @import foreach
#' @import iterators
#' @import data.table
#' @import grDevices
#' @import graphics
## usethis namespace: start
#' @importFrom lifecycle deprecate_soft
## usethis namespace: end
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate rename if_else select
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats model.frame model.matrix model.extract na.fail sd reshape
#' @description
#'
#' \lifecycle{maturing}
#'
#' Functions for memory-based learning
#' \if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}
#'
#' @details
#'
#' This is the version \code{2.0} (\code{'gordillo'}) of the package. It
#' implements a number of \code{R} functions useful for
#' modeling complex spectral spectra (e.g. NIR, IR).
#' The package includes functions for dimensionality reduction,
#' computing spectral dissimilarity matrices, nearest neighbor search,
#' and modeling spectral data using memory-based learning. 
#' 
#' Development versions can be found in the github repository of the package 
#' at \href{https://github.com/l-ramirez-lopez/resemble}{https://github.com/l-ramirez-lopez/resemble}.
#'
#' The functions available for dimensionality reduction are:
#' \itemize{
#'   \item{\code{\link{ortho_projection}}}
#'   \item{\code{\link{pc_projection}}}
#'   \item{\code{\link{pls_projection}}}
#'   \item{\code{\link{predict.ortho_projection}}}
#'   }
#' The functions available for computing dissimilarity matrices are:
#' \itemize{
#'   \item{\code{\link{dissimilarity}}}
#'   \item{\code{\link{f_diss}}}
#'   \item{\code{\link{cor_diss}}}
#'   \item{\code{\link{sid}}}
#'   \item{\code{\link{ortho_diss}}}
#'   }
#' The functions available for evaluating dissimilarity matrices are:
#' \itemize{
#'   \item{\code{\link{sim_eval}}}
#'   }
#' The functions available for nearest neighbor search:
#' \itemize{
#'   \item{\code{\link{search_neighbors}}}
#'   }
#' The functions available for modeling spectral data:
#' \itemize{
#'   \item{\code{\link{mbl}}}
#'   \item{\code{\link{mbl_control}}}
#'   }
#' Other supplementary functions:
#' \itemize{
#'   \item{\code{\link{plot.mbl}}}
#'   \item{\code{\link{plot.ortho_projection}}}
#'   }
#' @docType package
#' @name resemble-package
#' @aliases resemble-package resemble
#' @title Overview of the functions in the resemble package
#' @author 
#' 
#' \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} <ramirez.lopez.leo@@gmail.com> [aut, cre]
#' 
#' \href{https://orcid.org/0000-0002-1588-7519}{Antoine Stevens} [ctb] 
#' 
#' \href{https://orcid.org/0000-0003-1540-4748}{Raphael Viscarra Rossel} [ctb] 
#' 
#' \href{https://orcid.org/0000-0001-5416-8640}{Craig Lobsey} [ctb] 
#' 
#' \href{https://orcid.org/0000-0001-7325-9716}{Alex Wadoux} [ctb] 
#' 
#' \href{https://orcid.org/0000-0001-5695-8064}{Timo Breure} [ctb] 
######################################################################
# resemble
# Copyright (C) 2020 Leonardo Ramirez-Lopez
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
