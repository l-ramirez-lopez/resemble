#' @title Euclidean, Mahalanobis and cosine dissimilarity measurements
#' @description
#' This function is used to compute the dissimilarity between observations based on Euclidean or Mahalanobis distance measures or on cosine dissimilarity measures (a.k.a spectral angle mapper). 
#' @usage 
#' fDiss(Xr, X2 = NULL, method = 'euclid', center = TRUE, scaled = TRUE)
#' @param Xr a \code{matrix} (or \code{data.frame}) containing the (reference) data.
#' @param X2 an optional \code{matrix} (or \code{data.frame}) containing data of a second set of observations(samples).
#' @param method the method for computing the dissimilarity matrix. Options are \code{'euclid'} (Euclidean distance), \code{'mahalanobis'} (Mahalanobis distance) and \code{'cosine'} (cosine distance, a.k.a spectral angle mapper).
#' @param center a logical indicating if the spectral data \code{Xr} (and \code{X2} if specified) must be centered. If X2 is specified the data is scaled on the basis of \eqn{Xr \cup Xu}.
#' @param scaled a logical indicating if \code{Xr} (and \code{X2} if specified) must be scaled. If X2 is specified the data is scaled on the basis of \eqn{Xr \cup Xu}.
#' @details
#' In the case of both the Euclidean and Mahalanobis distances, the dissimilarity matrix \eqn{D} between between samples in a given matrix \eqn{X} is computed as follows:
#' \deqn{
#'      D(x_i, x_j) = \sqrt{(x_i - x_j)M^{-1}(x_i - x_j)^{\mathrm{T}}}
#'      }
#' where \eqn{M} is the identity matrix in the case of the Euclidean distance and the variance-covariance matrix of \eqn{M} in the case of the Mahalanobis distance.
#' The Mahalanobis distance can also be viewed as the Euclidean distance after applying a linear transformation of the original variables. 
#' Such a linear transformation is carried by using a factorization of the inverse covariance matrix as \eqn{M^{-1} = W^{\mathrm{T}}W}, where \eqn{M} is merely the square root of \eqn{M^{-1}} which can be found by using a singular value decomposition. 
#' Note that when attempting to compute the Mahalanobis distance on a dataset with highly correlated variables (i.e. spectral variables) the variance-covariance matrix may result in a singular matrix which cannot be inverted and therefore the distance cannot be computed. 
#' This is also the case when the number of samples in the dataset is smaller than the number of variables. 
#' For the computation of the Mahalanobis distance, the mentioned method is used. 
#' On the other hand the cosine dissimilarity \eqn{S} between two obsvervations \eqn{x_i} and \eqn{x_j} is computed as follows:
#' \deqn{
#' S(x_i, x_j) = cos^{-1}{\frac{\sum_{k=1}^{p}x_{i,k} x_{j,k}}{\sqrt{\sum_{k=1}^{p} x_{i,k}^{2}} \sqrt{\sum_{k=1}^{p} x_{j,k}^{2}}}}
#'      }
#' where \eqn{p} is the number of variables of the observations.
#' @return 
#' a \code{matrix} of the computed dissimilarities. 
#' @author Leonardo Ramirez-Lopez and Antine Stevens
#' @examples
#' require(prospectr)
#' 
#' data(NIRsoil)
#' 
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' # Euclidean distances between all the samples in Xr
#' ed <- fDiss(Xr = Xr, method = "euclid", center = TRUE, scaled = TRUE)
#' 
#' # Euclidean distances between samples in Xr and samples in Xu
#' ed.xr.xu <- fDiss(Xr = Xr, X2 = Xu, method = "euclid", center = TRUE, scaled = TRUE)
#' 
#' # Mahalanobis distance computed on the first 20 spectral variables
#' md.xr.xu <- fDiss(Xr = Xr[,1:20], X2 = Xu[,1:20], method = "mahalanobis", center = TRUE, scaled = TRUE)
#' 
#' # Cosine dissimilarity matrix
#' cdiss.xr.xu <- fDiss(Xr = Xr, X2 = Xu, method = "cosine", center = TRUE, scaled = TRUE)
#' @export

fDiss <- function(Xr, X2 = NULL, method = "euclid", center = TRUE, scaled = TRUE)
{
  if(!is.null(X2))
    if(ncol(X2) != ncol(Xr))
      stop("The number of columns (variables) in Xr must be equal to the number of columns (variables) in X2")

  mtd <- match.arg(method, c("euclid", "mahalanobis","cosine"))
  if(length(mtd) > 1)
    message(paste("More than one method was specified, only", mtd, "was used."))
  
  if(!is.logical(center))
    stop("'center' argument must be logical")
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical")
  
  if(center | scaled | method == "mahalanobis")
  {
    X <- rbind(Xr, X2)
    
    if(center){
      X <- sweep(x = X, MARGIN = 2, FUN = "-", STATS = colMeans(X))
    }
    
    if(scaled){
      X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = cSds(X))
    }
    
    if(method == "mahalanobis")
    {
      if(nrow(X) < ncol(X))
        stop("For computing the Mahalanobis distance, the total number of samples (rows) must be larger than the number of variables (columns).")
      X <- try(e2m(X, sm.method = "svd"), TRUE)
      if(!is.matrix(X))
        stop("The covariance matrix (for the computation of the Mahalanobis distance) is exactly singular. Try another metod.")
      method <- "euclid"
    }
      
    if(!is.null(X2))
    {
      X2 <- X[(nrow(X) - nrow(X2) +1):nrow(X), , drop = FALSE]
      Xr <- X[1:(nrow(X) - nrow(X2)), ]
    }else{
      Xr <- X
    }
    rm(X)
  }
  
  if(!is.null(X2))
  {
    dstc <- fastDist(X2, Xr, method)
    ptnc <- ifelse(method == "euclid", 0.5, 1)
    dvn <- ifelse(method == "euclid", ncol(Xr), 1)
    rslt <- (dstc^ptnc)/dvn
    colnames(rslt) <- paste("X2", 1:nrow(X2), sep = ".")
    rownames(rslt) <- paste("Xr", 1:nrow(Xr), sep = ".")
  }else{
    dstc <- fastDist(Xr, Xr, method)
    ptnc <- ifelse(method == "euclid", 0.5, 1)
    dvn <- ifelse(method == "euclid", ncol(Xr), 1)
    rslt <- (dstc^ptnc)/dvn
    rownames(rslt) <- paste("Xr", 1:nrow(Xr), sep = ".")
    colnames(rslt) <- rownames(rslt)
  }
  return(rslt)
}

#' @title A function for transforming a matrix from its Euclidean spcae to its Mahalanobis space
#' @keywords internal
e2m <- function(X, sm.method = c("svd", "eigen")){
  
  nms <- dimnames(X)
  
  if(ncol(X) > nrow(X))
    stop("In order to project the matrix to a Mahalanobis space, the number of observations of the input matrix must greater than its number of variables")
  
  if(length(sm.method) > 1)
    sm.method <- sm.method[1]
  if(!(sm.method %in% c("svd", "eigen")))
    stop("sm.method must be one of 'svd', 'eigen'")
  
  X <- as.matrix(X)
  vcv <- cov(X)
  sq_vcv <- sqrtSm(vcv, method = sm.method)
  sq_S <- solve(sq_vcv)
  ms_x <- X %*% sq_S
  dimnames(ms_x) <- nms
  return(ms_x)
}

#' @title Square root of (square) symetric matrices
#' @keywords internal
sqrtSm <- function(X, method = c("svd", "eigen")){
  
  if(!isSymmetric(X))
    stop("X must be a square symmetric matrix")
  if(length(method) > 1)
    method <- method[1]
  if(!(method %in% c("svd", "eigen")))
    stop("method must be one of 'svd', 'eigen'")
  
  if(method == "svd")
  {
    out <- svd(X)
    D <- diag(out$d)
    U <- out$v
    return(U %*% sqrt(D) %*% t(U))
  }
  
  if(method == "eigen")
  {
    out <- eigen(X)
    D <- diag(out$values)
    U <- out$vectors
    return(U %*% sqrt(D) %*% t(U))
  }
  
#   if(method == "Schur")
#   {
#     
#     out <- Schur(X)
#     D <- diag(out$EValues)
#     U <- out$Q
#     return(U %*% sqrt(D) %*% t(U))
#   }
}