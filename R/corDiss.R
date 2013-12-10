#' @title correlation and moving correlation dissimilarity measurements (corDiss)
#' @description
#' Computes correlation and moving correlation dissimilarity matrices. 
#' @usage 
#' corDiss(Xr, X2 = NULL, ws = NULL, center = TRUE, scaled = TRUE)
#' @param Xr a \code{matrix} (or \code{data.frame}) containing the (reference) data.
#' @param X2 an optional \code{matrix} (or \code{data.frame}) containing data of a second set of observations(samples).
#' @param ws for moving correlation dissimilarity, an odd integer value which specifies the window size. If \code{ws = NULL}, then the window size will be equal to the number of variables (columns), i.e. instead moving correlation, the normal correlation will be used. See details.
#' @param center a logical indicating if the spectral data \code{Xr} (and \code{X2} if specified) must be centered. If X2 is specified the data is scaled on the basis of \eqn{Xr \cup Xu}.
#' @param scaled a logical indicating if \code{Xr} (and \code{X2} if specified) must be scaled. If X2 is specified the data is scaled on the basis of \eqn{Xr \cup Xu}.
#' @details
#' The correlation dissimilarity \eqn{cd} between two obsvervations \eqn{x_i} and \eqn{x_j} is computed as follows:
#' \deqn{
#'      cd(x_i, x_j) = \frac{1}{2}(1 - cor(x_i, x_j))
#'      }
#' The avobe formlula is used when \code{ws = NULL}.
#' On the other hand (when \code{ws != NULL}) the moving correlation dissimilarity \eqn{mcd} between two obsvervations \eqn{x_i} and \eqn{x_j} is computed as follows:
#' \deqn{
#'      mcd(x_i, x_j) = \frac{1}{2 ws}\sum_{k=1}^{p-ws}(1 - cor(x_{i,(k:k+ws)}, x_{j,(k:k+ws)}))
#'      }
#' where \eqn{ws} represents a given window size which rolls sequantially fom 1 up to \eqn{p - ws} and  \eqn{p} is the number of variables of the observations.
#' @return 
#' a \code{matrix} of the computed dissimilarities. 
#' @author Antine Stevens and Leonardo Ramirez-Lopez
#' @examples
#' require(prospectr)
#' data(NIRsoil)
#' 
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' corDiss(Xr = Xr)
#' 
#' corDiss(Xr = Xr, X2 = Xu)
#'
#' corDiss(Xr = Xr, ws = 41)
#' 
#' corDiss(Xr = Xr, X2 = Xu, ws = 41) 
#' @export

corDiss <- function(Xr, X2 = NULL, ws = NULL, center = TRUE, scaled = TRUE)
{
  if(!is.null(X2))
    if(ncol(X2) != ncol(Xr))
      stop("The number of columns (variables) in Xr must be equal to the number of columns (variables) in X2")
  
  if(!is.logical(center))
    stop("'center' argument must be logical")
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical")
  
  if(center | scaled)
  {
    X <- rbind(Xr, X2)
    
    if(center){
      X <- sweep(x = X, MARGIN = 2, FUN = "-", STATS = colMeans(X))
    }
    
    if(scaled){
      X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = cSds(X))
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
  
  if(!is.null(ws))
  {
    if(ws < 3 | length(ws) != 1) 
      stop(paste("'ws' must be an unique odd value greater than 2")) 
    if((ws %% 2) == 0)
      stop("'ws' must be an odd value")
    if(ws >= ncol(Xr))
      stop("'ws' must lower than the number of columns (variables) in Xr")
    if(!is.null(X2))
    {
      rslt <- movcorDist(X2, Xr, ws)
      colnames(rslt) <- paste("X2", 1:nrow(X2), sep = ".")
      rownames(rslt) <- paste("Xr", 1:nrow(Xr), sep = ".")
    }else{
      rslt <- movcorDist(Xr, Xr, ws)
      colnames(rslt) <- paste("Xr", 1:nrow(Xr), sep = ".")
      colnames(rslt) <- rownames(rslt)
    }
  }else{
    if(!is.null(X2))
    {
      rslt <- fastDist(X2, Xr, "cor")
      colnames(rslt) <- paste("X2", 1:nrow(X2), sep = ".")
      rownames(rslt) <- paste("Xr", 1:nrow(Xr), sep = ".")
    }else{
      rslt <- fastDist(Xr, Xr, "cor")
      colnames(rslt) <- paste("Xr", 1:nrow(Xr), sep = ".")
      colnames(rslt) <- rownames(rslt)
    }
  }
  return(rslt)
}