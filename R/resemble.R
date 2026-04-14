#' @useDynLib resemble
#' @import lifecycle
#' @import Rcpp
#' @import foreach
#' @import iterators
#' @import grDevices
#' @import graphics
#' @import mathjaxr
## usethis namespace: start
#' @importFrom lifecycle deprecate_soft
## usethis namespace: end
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats model.frame model.matrix model.extract na.fail sd reshape
#' @importFrom stats approx median na.omit predict setNames
#' @importFrom utils flush.console tail
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads
#' @importFrom foreach foreach %dopar%
#' @description
#' \ifelse{html}{\out{<a href='https://www.tidyverse.org/lifecycle/#maturing'><img src='figures/lifecycle-maturing.svg' alt='Maturing lifecycle'></a>}}{\strong{Maturing}}
#'
#' Functions for memory-based learning and spectra search.
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
#' upon the methods presented in 
#' Ramirez-Lopez et al. (2013a) \doi{10.1016/j.geoderma.2012.12.014}, 
#' Ramirez-Lopez et al. (2026a) and Ramirez-Lopez et al. (2026b).
#'
#' Development versions can be found in the github repository of the package
#' at \href{https://github.com/l-ramirez-lopez/resemble}{https://github.com/l-ramirez-lopez/resemble}.
#'
#' The functions available for computing dissimilarity matrices are:
#' \itemize{
#'   \item{\code{\link{dissimilarity}}}: Computes a dissimilarity matrix based on a specified method.
#'   \item{\code{\link{diss_pca}}}: constructor for principal components-based dissimilarity method.
#'   \item{\code{\link{diss_pls}}}: constructor for partial least squares-based dissimilarity method.
#'   \item{\code{\link{diss_correlation}}}: constructor for correlation-based dissimilarity method.
#'   \item{\code{\link{diss_euclidean}}}: constructor for euclidean distance-based dissimilarity method.
#'   \item{\code{\link{diss_mahalanobis}}}: constructor for Mahalanobis distance-based dissimilarity method.
#'   \item{\code{\link{diss_cosine}}}: constructor for cosine-based dissimilarity method.
#'   }
#' The functions available for evaluating dissimilarity matrices are:
#' \itemize{
#'   \item{\code{\link{diss_evaluate}}}: Evaluates the effectiveness of a dissimilarity matrix using side information.
#'   }
#' The functions available for nearest neighbor search:
#' \itemize{
#'   \item{\code{\link{search_neighbors}}}: Search for nearest neighbors of a query spectrum in a reference dataset based on a specified dissimilarity method.
#'   }
#' The functions available for modeling spectral data:
#' \itemize{
#'   \item \code{\link{mbl}}: Memory-based learning for modeling spectral data.
#'   \item \code{\link{gesearch}}: An evolutionary method to search optimal samples in large spectral datasets.
#'   \item \code{\link{liblex}}: Builds a library of reusable localized models.
#'   }
#' The functions available for dimensionality reduction are:
#' \itemize{
#'   \item{\code{\link{ortho_projection}}}: Computes an orthogonal projection of the data based on either principal components or partial least squares.
#'   \item{\code{\link{predict.ortho_projection}}}
#'   }
#' @name resemble-package
#' @aliases resemble-package resemble
#' @title Overview of the functions in the resemble package
#' @references
#'
#' Ramirez-Lopez, L., Viscarra Rossel, R., Behrens, T., Orellano, C.,
#' Perez-Fernandez, E., Kooijman, L., Wadoux, A. M. J.-C., Breure, T.,
#' Summerauer, L., Safanelli, J. L., & Plans, M. (2026a). When spectral
#' libraries are too complex to search: Evolutionary subset selection for
#' domain-adaptive calibration. \emph{Analytica Chimica Acta}, under review.
#' 
#' Ramirez-Lopez, L., Metz, M., Lesnoff, M., Orellano, C.,
#' Perez-Fernandez, E., Plans, M., Breure, T., Behrens, T.,
#' Viscarra Rossel, R., & Peng, Y. (2026b). Rethinking local spectral
#' modelling: From per-query refitting to model libraries. 
#' \emph{Analytica Chimica Acta}, under review.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
#' Scholten, T. (2013a). The spectrum-based learner: A new local approach for
#' modeling soil vis-NIR spectra of complex data sets. \emph{Geoderma} 195-196,
#' 268-279. 
#'
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
#' \item Antoine Stevens (\href{https://orcid.org/0000-0002-1588-7519}{ORCID})
#' \item Claudio Orellano
#' }
######################################################################
# resemble
# Copyright (C) 2014-2026 Leonardo Ramirez-Lopez
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
"_PACKAGE"
NULL
