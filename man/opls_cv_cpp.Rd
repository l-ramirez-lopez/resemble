% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{opls_cv_cpp}
\alias{opls_cv_cpp}
\title{Internal Cpp function for performing leave-group-out cross-validations for pls regression}
\usage{
opls_cv_cpp(X, Y, scale, method, 
                  mindices, pindices, 
                  min_component, ncomp, 
                  new_x, 
                  maxiter, tol, 
                  wapls_grid)
}
\arguments{
\item{X}{a matrix of predictor variables.}

\item{Y}{a matrix of a single response variable.}

\item{scale}{a logical indicating whether the matrix of predictors (\code{X}) must be scaled.}

\item{method}{the method used for regression. One of the following options: \code{'pls'} or \code{'wapls'} or \code{'completewapls1p'}.}

\item{mindices}{a matrix with \code{n} rows and \code{m} columns where \code{m} is equivalent to the number of 
resampling iterations. The elements of each column indicate the indices of the observations to be used for modeling at each 
iteration.}

\item{pindices}{a matrix with \code{k} rows and \code{m} columns where \code{m} is equivalent to the number of 
resampling iterations. The elements of each column indicate the indices of the observations to be used for predicting at each 
iteration.}

\item{min_component}{an integer indicating the number of minimum pls components (if the \code{method = 'pls'}).}

\item{ncomp}{an integer indicating the number of pls components.}

\item{new_x}{a matrix of one row corresponding to the observation to be predicted (if the \code{method = 'wapls'}).}

\item{maxiter}{maximum number of iterations.}

\item{tol}{limit for convergence of the algorithm in the nipals algorithm.}

\item{wapls_grid}{the grid on which the search for the best combination of minimum and maximum pls factors of \code{'wapls'} is based on in case \code{method = 'completewapls1p'}.}
}
\value{
a list containing the following one-row matrices:
\itemize{
\item{\code{rmse_seg}}{ the RMSEs.}
\item{\code{st_rmse_seg}}{ the standardized RMSEs.}
\item{\code{rsq_seg}}{ the coefficients of determination.}
}
}
\description{
For internal use only!.
}
\author{
Leonardo Ramirez-Lopez
}
\keyword{internal}
