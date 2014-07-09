#-----------------------------------------------
# Load packages and prepare data
#-----------------------------------------------

# load packages
# Don't forget to install required packages if needed...
library(resemble);library(prospectr)
library(ggplot2)
library(foreach); library(reshape2)

# load dataset from the "Chimiométrie 2006" challenge
# see ?NIRsoil for an explanation of the variables
data(NIRsoil) 

# Filter the data using the Savitzky and Golay smoothing filter with a window size of
# 11 spectral variables and a polynomial order of 3 (no differentiation).
NIRsoil$spc_sg <- savitzkyGolay(X = NIRsoil$spc, p = 3, w = 11, m = 0)

# bin spectra every 10 nm 
NIRsoil$spc_sg <- binning(X=NIRsoil$spc, bin.size=5)
wavelength <- as.numeric(colnames(NIRsoil$spc)) # store wavelength position

# Select dependent variable
dep <- "Nt" # Total Nitrogen

# Remove NA values in dependent variable
anyna <- is.na(NIRsoil[,dep])
NIRsoil <- NIRsoil[!anyna,]

# Select cal & val dataset
Yr <- NIRsoil[as.logical(NIRsoil$train),dep] # dependent variable (calibration set)
Xr <- NIRsoil$spc_sg[as.logical(NIRsoil$train),] # spectral matrix (calibration set)
Yu <- NIRsoil[!as.logical(NIRsoil$train),dep] # dependent variable (validation set)
Xu <- NIRsoil$spc_sg[!as.logical(NIRsoil$train),] # spectra matrix (validation set)

#-----------------------------------------------
# Computing spectral dissimilarities
#-----------------------------------------------

# Compute distances
euDis <- fDiss(Xr = Xr,method="euclid") # euclidean distance
cosineDis <- fDiss(Xr = Xr, method="cosine") # cosine distance
corDis <- corDiss(Xr = Xr) # correlation distance
mcorDis <- corDiss(Xr = Xr, ws = 41) # moving correlation distance
pcaDis <- orthoDiss(Xr,pcSelection=list("cumvar",.99),method="pca",local=F)$dissimilarity # default parameter in orthoDiss
# Computation of a principal component dissimilarity matrix using the
# "opc" method for the selection of the principal components
opcDis <- orthoDiss(Xr,Yr = Yr,pcSelection = list("opc", 40),method = "pca")$dissimilarity
# Computation of a partial least squares (PLS) dissimilarity matrix 
plsDis <- orthoDiss(Xr,Yr = Yr,pcSelection = list("manual", 10),method = "pls")$dissimilarity

# Select one sample with the smallest mahalanobis distance to the centre of the data
# First Project in the PC space
pca <- prcomp(Xr,.center=T) # pca...
sc <- scale(pca$x,center=T,scale=T) # standardized score
pvar <- cumsum(pca$sdev^2/sum(pca$sdev^2)) # cumulative explained variance
pc <- max(which(pvar < .99)) + 1 # number of PC's accounting at least for .99 % of the variation
sc <- sc[,1:pc] # truncated score matrix
mahal <- mahalanobis(sc,center=colMeans(sc),cov=cov(sc)) # mahalanobis distance
id <- which.min(mahal) # index of the sample closest to the centre of the data

# Plot the closest sample to the mean
p <- ggplot(data=melt(Xr),aes(x=as.numeric(as.character(Var2)),y=value,group=Var1)) + geom_line(col="grey80",alpha=.5) + theme_bw() + labs(x="wavelength /nm",y="Absorbance") + scale_x_continuous(breaks=seq(1000,2500,250))
p + geom_line(data=melt(Xr[id,,drop=F]),col="red")

# Get distances to the selected sample and stack data
d <- as.data.frame(cbind(euDis[,id],cosineDis[,id],corDis[,id],mcorDis[,id],pcaDis[,id],opcDis[,id],plsDis[,id]))
colnames(d) <- paste0(c("euclidean","cosine","correlation","moving correlation","mahalanobis","opc","opc with pls")," distance")

# for each dissimilarity measure, find the 10 closest neighbours to the selected sample
knn <- sapply(d,function(x)order(x))[1:11,]
knn

# Plot
# merge spectra into one data.frame and format for ggplot2
spc <- foreach(i = 1:ncol(knn),.combine=rbind)%do%{
  data.frame(melt(Xr[knn[,i],]),d=colnames(knn)[i])
}
p <- ggplot(data=spc,aes(x=as.numeric(as.character(Var2)),y=value,group=Var1)) + geom_line(col="grey80") + facet_wrap(~d) + theme_bw() + labs(x="wavelength /nm",y="Absorbance") + scale_x_continuous(breaks=seq(1000,2500,250))
p + geom_line(data=melt(Xr[id,,drop=F]),col="red")

# Now we evaluate whether dissimarities are representative of the sample compositional variation
# using simEval
euSE <- simEval(d = euDis, sideInf = Yr)
cosineSE <- simEval(d = cosineDis, sideInf = Yr)
corSE <- simEval(d = corDis, sideInf = Yr)
mcorSE <- simEval(d = mcorDis, sideInf = Yr)
pcaSE <- simEval(d = pcaDis, sideInf = Yr)
opcSE <- simEval(d = opcDis, sideInf = Yr)
plsSE <- simEval(d = plsDis, sideInf = Yr)
sim <- list(euSE,cosineSE,corSE,mcorSE,pcaSE,opcSE,plsSE)
names(sim) <-  paste0(c("euclidean","cosine","correlation","moving correlation","mahalanobis","opc","opc with pls")," distance")

# stack rmsd and r (correlation) and plot
simplot <- (do.call(rbind,lapply(sim,function(x)x$eval)))
simplot$d <- factor(rownames(simplot),levels=rownames(simplot))
ggplot(data=melt(simplot),aes(y=value,x=d,fill=d)) + geom_bar(stat="identity") + facet_wrap(~variable,scale="free_y") + theme_bw() + labs(x="dissimilarities",y="") + scale_fill_brewer("Dissimilarities",palette="Set1") + scale_x_discrete(labels=abbreviate) + theme(legend.position="top")
# the 'opc' distance is the most representative of the soil compositional variation
# this is expected since 'opc' minimizes the distance between the closest neighbours and the training samples in the compositional space

#-----------------------------------------------
# Predict sample compostion with 'mbl'
# see ?mbl for other examples
#-----------------------------------------------

# Ex 1 : A simple Memory-based learning approach using partial least square regression and correlation distance
# First, set the parameters controlling the mbl function
ctrl <- mblControl(sm = "cor",              # dissimilarity
                   center = TRUE,           # Does the predictors need to be centered ?
                   scaled = TRUE,           # Does the predictors need to be scaled ?
                   valMethod = "none",      # No internal validation method
                   range.pred.lim = TRUE,   # Are the predictions constrained by the range of the response variable in each local model ?
                   progress = TRUE          # progress bar
                   )   
# Run mbl
mbl_cor <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl,
               dissUsage = "none",   # How the distance are used: 'none' = only in the selection of neighbours, 'weights' = distances used to computed sample weight through the tricubic function
               k = seq(25, 150, by = 25),  # a vector of number of neighbours to be tested
               method = "wapls1", # regression method: the prediction is a weighted average of values predicted with increasing pls components
               pls.c = c(4,15) # number of pls components                 
               )
mbl_cor
plot(mbl_cor) # plot of the error vs k

# Get predictions and plot pred-obs
preds <- getPredictions(mbl_cor)
predobs <- data.frame(pred=melt(preds),obs=Yu) 
predobs$k <- sub(".+_([0-9]+)","\\1",predobs$pred.variable)
p <- ggplot(data=predobs,aes(x=pred.value,y=obs))  + geom_point() + facet_wrap(~k) + theme_bw() + geom_abline(col="red") + labs(x="Predicted Nt",y="Observed Nt")
p 
p + geom_text(data=mbl_cor$YuPredictionStats,aes(x=-Inf,y=Inf,label=paste0("rmse = ",round(rmse,2))),vjust=2,hjust=-0.1) # this adds rmse info

# Ex 2 : same as  Ex. 1 but with the pc distance
ctrl$sm  <- "pc"
ctrl$pcSelection <- list(method="cumvar",value=.999)  # select as many PC's as to explain 99.9% of the spectral variation

mbl_pc <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl,
               dissUsage = "none",   
               k = seq(25, 150, by = 25), 
               method = "wapls1", 
               pls.c = c(4,15)
              )
mbl_pc

# Ex 3 : same as Ex. 2 but with the pls distance
ctrl$sm  <- "pls"
ctrl$pcSelection <- list(method="manual",value=10)

mbl_pls <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl,
              dissUsage = "none",   
              k = seq(25, 150, by = 25),  
              method = "wapls1", 
              pls.c = c(4,15)
              )
mbl_pls # this shows some improvement compared to pc distance...

# Ex 4 : same as Ex. 2 but with a user-defined distance matrix
# distances could be for instance computed on the 1st derivative spectra
Xr_der1 <- savitzkyGolay(NIRsoil$spc[as.logical(NIRsoil$train),], p = 5, w = 21, m = 1)
Xu_der1 <- savitzkyGolay(NIRsoil$spc[!as.logical(NIRsoil$train),], p = 5, w = 21, m = 1)
pls_dis <- orthoDiss(Xr=Xr_der1,X2 = Xu_der1, Yr = Yr,method = "pls", pcSelection = list("manual",10)) # pls distance in the 1st derivative space

ctrl$sm  <- "none" # we use pls_dis instead...
mbl_user_diss <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,mblCtrl = ctrl,
                  dissimilarityM = pls_dis$dissimilarity, # user-specified dissimilarity matrix
                  dissUsage = "none",   
                  k = seq(50, 150, by = 25),  
                  method = "wapls1",  # we use the 'standard' pls algorithm...
                  pls.c = c(4,15)       # with a max of 15 components
                  )
mbl_user_diss

# Ex 5 : Let's try with another multivariate regression technique: pls
# the number of pls component for each local model is computed by cross-validation
ctrl$sm  <- "pls" 
mbl_pls_cv <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl,
               dissUsage = "none",   
               k = seq(25, 150, by = 25),  
               method = "pls",  # we use the 'standard' pls algorithm...
               pls.c = 15       # with a max of 15 components
               )
mbl_pls_cv # this shows some improvement compared to wapls1 regression

# Ex 6 : Let's try with another multivariate regression technique: gpr
mbl_gpr <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl,
               dissUsage = "none",   
               k = seq(25, 150, by = 25),  
               method = "gpr"  # predictions are done through gaussian process                   
                )
mbl_gpr # this is comparable with pls, without the need for tuning the number of components...

# Ex 7: A Memory-based learning approach (the spectrum-based learner), as implemeted in Ramirez-Lopez et al. (2013) 
sbl <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
           mblCtrl =  mblControl(valMethod = "none"), # most default values in mblControl correspond to the spectrum-based learner..., 
           dissUsage = "predictors", # the distance matrix is used as additional predictor
           k = seq(25, 150, by = 25), 
           method = "gpr" # predictions are done through gaussian process 
           )
sbl

#Plot all the results together
mymodels <- list(cor = mbl_cor,pc = mbl_pc, pls = mbl_pls, user_def = mbl_user_diss , pls_cv = mbl_pls_cv, gpr = mbl_gpr,sbl = sbl) 
tmp <- do.call(rbind,lapply(names(mymodels), function(x) data.frame(name = x, mymodels[[x]]$YuPredictionStats))) # stack model results
ggplot(data=tmp,aes(x=k,y=rmse,colour=name,linetype=name)) + geom_point() + geom_line() + theme_bw() +  ylim(c(0.4,0.8)) + scale_x_continuous(breaks=seq(25, 150, by = 25))

# Ex 7 : 
# In the previous examples, the selected samples for the local models (ie neighbours) are selected with a fixed number (k)
# but the neighborhood could be defined by distance tresholds
# this is implemented in mbl with the k.diss argument
mbl_kdiss <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, 
               mblCtrl = mblControl(sm = "cor",valMethod = "none"), # correlation distance
               dissUsage = "none",   
               k.diss = c(0.25,0.1,0.05),  # samples with a correlation less than 0.75, 0.9, 0.95 with the sample to predict will not be included in the models
               k.range = c(20,nrow(Xr)), # the minimum and maximum number of samples to be included in the local calibrations
               method = "gpr"  
              )
# let's see what is the number of samples used in the calibrating models
res <-mbl_kdiss$results[[2]] # prediction table for correlation distance = 0.1
hist(as.numeric(as.factor(res$k)),xlab="k",main="")

# setting the tresholds might be less intuitive for other distances
# we suggest to pre-compute the dissimilarity matrix
# and look at the histogram of the dissimilarity values
# here is an example with opc distance
opc_dis <- orthoDiss(Xr=Xr,X2 = Xu, Yr = Yr,method = "pca", pcSelection = list("opc",40)) # pls distance in the 1st derivative space
ctrl$sm  <- "none" 
mbl_kdiss <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl,
                 dissimilarityM = opc_dis$dissimilarity,                  
                 dissUsage = "none",   
                 k.diss = c(0.05,0.1,0.2,0.5),  # samples with a pc distance > 0.1, 0.2, 0.5  with the sample to predict will not be included in the models
                 k.range = c(20,nrow(Xr)), # the minimum and maximum number of samples to be included in the local calibrations
                 method = "gpr"  
                )

# Ex 8: Internal validation
# resemble provides a way to make an internal validation of the results: 
# 'NNv': the nearest neighbour to the sample to predict is left out of the list of k nearest neighbours and predicted and compared with the actual value to compute a rmse
# 'loc_crossval' : each local model is cross-validated and the mean of the cross-validation results is computed
ctrl <- mblControl(sm="pls",pcSelection = list("manual",10),
                   valMethod = c("loc_crossval","NNv"),
                   resampling = 10, # controls the local cross-validation: this is the number of partition,
                   p = 0.25 # percentage of samples included in each partition
                  )
mbl_valid <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl,
               dissUsage = "none",   
               k = seq(25, 150, by = 25),  
               method = "gpr"  
              ) # note that the algorithm is slower when using internal validation

plot(mbl_valid) # shows the external (against Yr) and internal rmse, which seems a bit optimistic, compared to external validation

# Note that to be correct, external validation should be done with 3 sets:
# a reference (Xr), tuning and test sets
# the tuning set is used to find the best combination of MBL parameters (k, sm, etc..)
# which is then applied to predict the test set

#-----------------------------------------------
# Example of parallel execution
#-----------------------------------------------
# To speed-up mbl, the local models can be run in parallel
# using the 'parallel' package and the parallel back-end for the foreach %dopar% function 
# provided by 'doParallel'
library(parallel)
library(doParallel)
# here is the time spent with one core
ctrl <- mblControl(sm="pls",pcSelection = list("manual",10),progress = T)
<<<<<<< HEAD
system.time(mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl, dissUsage = "none", k = seq(25, 150, by = 25),  method = "gpr"))
# now with 4 cores
cl <- makeCluster(4) # create a set 4 running R copies
registerDoParallel(cl) # register them to work with foreach
system.time(mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl, dissUsage = "none", k = seq(25, 150, by = 25),  method = "gpr"))
=======
system.time(mbl_test <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl, dissUsage = "none", k = seq(25, 150, by = 25),  method = "gpr"))
# now with 4 cores
cl <- makeCluster(4) # create a set 4 running R copies
registerDoParallel(cl) # register them to work with foreach
system.time(mbl_test <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl, dissUsage = "none", k = seq(25, 150, by = 25),  method = "gpr"))
>>>>>>> d9e8596248ed123ae8764316c71d1c81a06d911b
# I get a ~3 x speed-up ...
registerDoSEQ() # un-register
stopCluster(cl) # delete R instances