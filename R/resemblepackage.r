#'@description 
#' This package implements a number of \code{R} functions useful for modelling 
#' complex visible and infrared spectra(\acronym{vis-IR}). The packages includes functions for 
#' for projecting spectral data onto orthogonal spaces, computing spectral similarity/dissimilarity
#' matrices, removing irrelevant spectra from a reference set, and modelling spectral data
#' using memory-based learning.
#' 
#' The functions available for projecting the spectra are:
#' \itemize{
#'   \item{\code{\link{orthoProjection}}} 
#'   }
#' The functions available for computing similarity/dissimilarity matrices are:
#' \itemize{
#'   \item{\code{\link{fDiss}}} 
#'   \item{\code{\link{corDiss}}} 
#'   \item{\code{\link{sid}}} 
#'   \item{\code{\link{orthoDiss}}} 
#'   }
#' The functions available for evaluating similarity/dissimilarity matrices are:
#' \itemize{
#'   \item{\code{\link{simEval}}} 
#'   }
#' The functions available for removing irrelevant spectra from a reference set are:
#' \itemize{
#'   \item{\code{\link{neigCleaning}}} 
#'   }
#' The functions available for modelling spectral data using memory-based learning are:
#' \itemize{
#'   \item{\code{\link{mblController}}}
#'   \item{\code{\link{mbl}}} 
#'   }
#' Other supplementary functions are:
#' \itemize{
#'   \item{\code{\link{print.localOrthoDiss}}}
#'   \item{\code{\link{print.mbl}}}
#'   \item{\code{\link{plot.mbl}}}
#'   }
#'@docType package
#'@name resemble
#'@title resemble package
#'@import Rcpp RcppArmadillo foreach iterators 
#'@useDynLib resemble
#'@author Leonardo Ramirez-Lopez & Antoine Stevens
#'
NULL