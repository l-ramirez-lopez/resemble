#' @useDynLib resemble
#' @import lifecycle
#' @import Rcpp
#' @import foreach
#' @import iterators
#' @import data.table
#' @import grDevices
#' @import graphics
#' @import mathjaxr
## usethis namespace: start
#' @importFrom lifecycle deprecate_soft
## usethis namespace: end
#' @importFrom magrittr %>%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats model.frame model.matrix model.extract na.fail sd reshape
#' @description
#' \ifelse{html}{\out{<a href='https://www.tidyverse.org/lifecycle/#maturing'><img src='figures/lifecycle-maturing.svg' alt='Maturing lifecycle'></a>}}{\strong{Maturing}}
#'
#' Functions for memory-based learning
#'
#' \if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}
#'
#' @details
#' This is the version
#' `r paste(pkg_info()[1:2], collapse = " -- ")`
#' of the package. It implements a number of functions useful for
#' modeling complex spectral spectra (e.g. NIR, IR).
#' The package includes functions for dimensionality reduction,
#' computing spectral dissimilarity matrices, nearest neighbor search,
#' and modeling spectral data using memory-based learning. This package builds
#' upon the methods presented in Ramirez-Lopez et al. (2013) <doi:10.1016/j.geoderma.2012.12.014>>.
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
#' @references
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. 2013a. The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196,
#' 268-279.
#' @seealso
#' Useful links:
#' \itemize{
#' \item \url{https://github.com/l-ramirez-lopez/resemble}
#' \item Report bugs at \url{https://github.com/l-ramirez-lopez/resemble/issues}
#' }
#' @author
#'
#' \strong{Maintainer / Creator}: Leonardo Ramirez-Lopez \email{ramirez.lopez.leo@gmail.com}
#'
#' Authors:
#' \itemize{
#' \item Leonardo Ramirez-Lopez (\href{https://orcid.org/0000-0002-5369-5120}{ORCID})
#'
#' \item Antoine Stevens (\href{https://orcid.org/0000-0002-1588-7519}{ORCID})
#'
#' \item Claudio Orellano
#'
#' \item Raphael Viscarra Rossel (\href{https://orcid.org/0000-0003-1540-4748}{ORCID})
#'
#' \item Zefang Shen
#'
#' \item Craig Lobsey (\href{https://orcid.org/0000-0001-5416-8640}{ORCID})
#'
#' \item Alex Wadoux (\href{https://orcid.org/0000-0001-7325-9716}{ORCID})
#'
#' }
######################################################################
# resemble
# Copyright (C) 2022 Leonardo Ramirez-Lopez
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
