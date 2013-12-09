#' @title Plot method for an object of class \code{mbl}
#' @aliases plot.mbl
#' @usage \method{plot}{mbl}(x, g = c('validation', 'pca'), param = 'rmse', pcs = c(1,2), ...)
#' @param x an object of class \code{mbl} (as returned by \code{mbl}). 
#' @param g a character \code{vector} indicating what results shall be plotted. Options are: 'validation' (for plotting the validation results) and/or 'pca' (for plotting the principal components).
#' @param param one of the following options 'rmse', 'st.rmse' or 'r2'. The respective validation statistic is then plotted. It is only available if \code{'validation'} is specified in the \code{g} argument. 
#' @param pcs a vector of length one or two indicating the principal components to be plotted. Default is \code{c(1, 2)}. It is only available if \code{'pca'} is specified in the \code{g} argument. 
#' @param ... arguments to be passed to methods (not yet functional).
#' @authors Leonardo Ramirez-Lopez and Antoine Stevens
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematt?, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{mbl}}
#' @examples
#' \dontrun{
#' require(prospectr)
#' data(NIRsoil)
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' Xu <- Xu[!is.na(Yu),]
#' Yu <- Yu[!is.na(Yu)]
#' 
#' Xr <- Xr[!is.na(Yr),]
#' Yr <- Yr[!is.na(Yr)] 
#'
#' ctrl <- mblController(sm = "cor", ws = 51, pcSelection = list("cummvar", 0.999), valMethod = c("NNv"), scaled = TRUE, center = TRUE)
#'
#' ex1 <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
#'            mblCtrl = ctrl,
#'            distUsage = "none", 
#'            k = seq(30, 250, 30), 
#'            method = "wapls1",
#'            plsF = c(7, 20))
#'
#' plot(ex1)
#' }
#' @export

plot.mbl <- function(x, 
                     g = c("validation", "pca"), 
                     param = "rmse", 
                     pcs = c(1,2), ...){
  
  object <- x
  pm <- par()$mfrow

  if("validation" %in% g)
  {
    col <- NULL
    if(!is.null(object$nnValStats)){
      nnValStats <- cbind(object$nnValStats, val = "NNv")
      col <- c(col, "dodgerblue")
    } else {nnValStats <- NULL}
    if(!is.null(object$localCrossValStats)){
      localCrossValStats <- cbind(object$localCrossValStats, r2 = NA, val = "loc_crossval")
      col <- c(col, "green4")
    } else {localCrossValStats <- NULL}
    if(!is.null(object$YuPredictionStats)){
      YuPredictionStats <- cbind(object$YuPredictionStats, val = "Yu prediction")
      col <- c(col, "red")
    } else {YuPredictionStats <- NULL}
    
    tpl <- rbind(nnValStats, localCrossValStats, YuPredictionStats)
    
    if(is.null(tpl)){
      message("No validation results to plot")
    }else{
      par(mfrow = c(1, length(g)))
      dtn <- colnames(tpl) 
      opt <- c("rmse", "st.rmse", "r2")
      dt <- !is.element(dtn, opt[!is.element(opt, param)])
      toPlot <- reshape(tpl[,dt], timevar = "val", idvar = "k", direction = "wide")
      
      if(param == "r2"){
        toPlot <- toPlot[,!colnames(toPlot) == "r2.loc_crossval"]
        col <- col[!col == "green4"]
      }
      matplot(toPlot[,1], toPlot[,-1], 
              type="b", xlab = dtn[1], 
              ylab = param, pch = 1, ylim = c(min(toPlot[,-1]), 1.1 * max(toPlot[,-1])),
              col = col)
      mtext("Validation", col = "red")
      # Adding a legend
      legend("topright", legend = colnames(toPlot[,-1,drop = FALSE]), pch =1,
             col = col, cex = 0.8, border = "red",box.lty = 3, box.col = "grey")
      #is.element("ggplot2", installed.packages()[,"Package"])
    }
  }
  
  if("pca" %in% g)
  {
    if(object$pcAnalysis$n.componentsUsed == 1)
    {
      rng <- range(object$pcAnalysis$scores_Xr[,pcs[1]], object$pcAnalysis$scores_Xu[,pcs[1]])
      rng <- 1.2 * c(-max(abs(rng)), max(abs(rng)))
      tp <- c(object$pcAnalysis$scores_Xr[,1], object$pcAnalysis$scores_Xu[,1])
      tp <- data.frame(index = 1:length(tp), tp = tp, set = c(rep("Xr", length(object$pcAnalysis$scores_Xr[,1])), rep("Xu", length(object$pcAnalysis$scores_Xu[,1]))))
      tp <- tp[order(tp$tp),]
      tp$index <- 1: length(tp$index)
      plot(tp[tp$set == "Xr",1:2], ylim = rng, 
           col = rainbow(1, s = 1, v = 0, alpha = 0.3), 
           pch = 16, ylab = "pc1 (standardized)", xlab = "index (ordered values)")
      mtext("pc Analyisis", col = "red")
      points(tp[tp$set == "Xu",1:2], xlim = rng, ylim = rng, col = heat.colors(1, alpha = 0.6), pch = 16)
      legend("topright", legend = c("Xr", "Xu"),
             col = c(rainbow(1, s = 1, v = 0, alpha = 0.3), heat.colors(1, alpha = 0.4)), pch = 16, cex = 0.8, box.lty = 3, box.col = "grey")
    }else{
      rng <- range(object$pcAnalysis$scores_Xr[,pcs], object$pcAnalysis$scores_Xu[,pcs])
      rng <- 1.2 * c(-max(abs(rng)), max(abs(rng)))
      
      xl <- paste(colnames(object$pcAnalysis$scores_Xr[,pcs[1],drop=F]), " (standardized)", sep = "") 
      yl <- paste(colnames(object$pcAnalysis$scores_Xr[,pcs[2],drop=F]), " (standardized)", sep = "") 
       
      plot(object$pcAnalysis$scores_Xr[,pcs], xlab = xl, ylab = yl, xlim = rng, ylim = rng, 
           col = rainbow(1, s = 1, v = 0, alpha = 0.3), pch = 16)
      mtext("pc Analyisis", col = "red")
      
      points(object$pcAnalysis$scores_Xu[,pcs], xlim = rng, ylim = rng, col = heat.colors(1, alpha = 0.4), pch = 16)
      legend("topright", legend = c("Xr", "Xu"),
             col = c(rainbow(1, s = 1, v = 0, alpha = 0.3), heat.colors(1, alpha = 0.4)), pch = 16, cex = 0.8, box.lty = 3, box.col = "grey")
    }
  }
  par(mfrow = pm)
  title(main = "Memory-based learning results")
}