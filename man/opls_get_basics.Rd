% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{opls_get_basics}
\alias{opls_get_basics}
\title{fast orthogonal scores algorithn of partial leat squares (opls)}
\usage{
opls_get_basics(X, Y, ncomp, scale, 
                maxiter, tol, 
                algorithm = "pls", 
                xls_min_w = 3, 
                xls_max_w = 15)
}
\arguments{
\item{X}{a matrix of predictor variables.}

\item{Y}{a matrix of either a single or multiple response variables.}

\item{ncomp}{the number of pls components.}

\item{scale}{logical indicating whether \code{X} must be scaled.}

\item{maxiter}{maximum number of iterations.}

\item{tol}{limit for convergence of the algorithm in the nipals algorithm.}

\item{algorithm}{(for weights computation) a character string indicating
what method to use. Options are:
\code{'pls'} for pls (using covariance between X and Y),
\code{'mpls'} for modified pls (using correlation between X and Y) or
\code{'xls'} for extended pls (as implemented in BUCHI NIRWise PLUS software).}

\item{xls_min_w}{(for weights computation) an integer indicating the minimum window size for the "xls"
method. Only used if \code{algorithm = 'xls'}. Default is 3 (as in BUCHI NIRWise PLUS software).}

\item{xls_max_w}{(for weights computation) an integer indicating the maximum window size for the "xls"
method. Only used if \code{algorithm = 'xls'}. Default is 15 (as in BUCHI NIRWise PLUS software).}
}
\value{
a list containing the following elements:
\itemize{
\item{\code{coefficients}}{ the matrix of regression coefficients.}
\item{\code{bo}}{ a matrix of one row containing the intercepts for each component.}
\item{\code{Y_loadings}}{ the matrix of Y loadings.}
\item{\code{projection_mat}}{ the projection matrix.}
\item{\code{transf}}{ a \code{list} conating two objects: \code{Xcenter} and \code{Xscale}}.
}
}
\description{
Computes orthogonal socres partial least squares (opls)
regressions with the NIPALS algorithm. It allows multiple response variables.
In contrast to code{opls} function, this one does not compute unnecessary
data for (local) regression.
For internal use only!
}
\author{
Leonardo Ramirez-Lopez
}
\keyword{internal}
