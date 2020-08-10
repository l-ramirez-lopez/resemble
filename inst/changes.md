# Changes implemented for version 2.0


In general, the code has been reformatted for compliying with Google and 
tidyverse style guideleines. Here the major changes are listed.

### Main new fucntions
- search_neighbors
- dissimilarity


### orthoProjection, pcProjection, plsProjection (renamed to ortho_projection, 
pc_projection, pls_projection respectively):
- X2 argument renamed to Xu (for consistency throughout all the fucntions)
- Argument scaled renamed to .scale
- Argument max.iter renamed to max_iter
- Bug fix: when the "pca.nipals method was used and the method to select the pcs wa "opc",
  the function was returning 1 component less than the maximum requested.
- "pca.nipals" is now implemented in C++ via Rcpp
- Bug fix in plsProjection: when "cumvar" was used as the pcSelection method, an 
  error about data allocation in a matrix was retrieved
- Argument pcSelection renamed to pc_selection
- ... is deprecated in both pcProjection and plsProjection
- Argument cores is deprecated, it wwas used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. Multi-threading via OpenMP is still 
available, but the number of threads must be set in the R global enviroment (e.g.
Sys.setenv("OMP_NUM_THREADS" = <a integer>)).
- Names the following outputs have been changed: X.loadings, Y.loadings, sc.sdv 
and n.components. Their new names are: X_loadings, Y_loadings, sc_sdv and 
n_components.


### corDiss (renamed to cor_diss):
- X2 argument renamed to Xu (for consistency throughout all the fucntions)
- Argument scaled renamed to .scale
- default for .scale has changed from TRUE to FALSE
- the dimnames of the resulting matrix are now Xr_1... Xr_n (previusly Xr.1... Xr.n)

### fDiss (renamed to f_diss):
- X2 argument renamed to Xu (for consistency throughout all the fucntions)
- Argument scaled renamed to .scale
- default for .scale has changed from TRUE to FALSE
- the dimnames of the resulting matrix are now Xr_1... Xr_n (previusly Xr.1... Xr.n)
- argument method changed to diss_method

### sid:
- X2 argument renamed to Xu (for consistency throughout all the fucntions)
- Argument scaled renamed to .scale
- default for .scale has changed from TRUE to FALSE
- the dimnames of the resulting matrix are now Xr_1... Xr_n (previusly Xr.1... Xr.n)


### orthoDiss (renamed to ortho_diss):
- X2 argument renamed to Xu (for consistency throughout all the fucntions)
- Argument scaled renamed to .scale
- Argument local renamed to .local
- Argument pcSelection renamed to pc_selection
- Argument return.all renamed to compute_all
- Argument cores is deprecated, it wwas used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. Multi-threading via OpenMP is still 
available, but the number of threads must be set in the R global enviroment (e.g.
Sys.setenv("OMP_NUM_THREADS" = <a integer>)).
- When \code{.local = TRUE} a new output is produced: 'neighborhood_info' which 
is a data.frame containing the relevant information about the neighborhood of 
each sample (e.g. neighborhood indices, number of components used at each 
neighborhood, etc)
- Output global.variance.info has been renamed to global_variance_info


### simEval (renamed to sim_eval):
- argument sideInf renamed to side_info
- argument lower.tri renamed to lower_triangle
- argument cores renamed to omp_threads
- lower_triangle is deprecated. Now if a vector is passed to d, the function assumes
  that it is a lower triangle of a distance matrix
- Argument cores is deprecated, it wwas used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. Multi-threading via OpenMP is still 
available, but the number of threads must be set in the R global enviroment (e.g.
Sys.setenv("OMP_NUM_THREADS" = <a integer>)).
  
  
### mbl
- pls.max.iter, pls.tol and noise.v were moved to mbl from mbl_control()
- Argument scaled (from mbl_control()) renamed to .scale and moved to mbl
- new arguments: gh and spike
- Argument pcSelection renamed to pc_selection
- Argument mblCtrl renamed to control
- Argument dissUsage renamed to diss_usage
- order of the Yr, Xr, Yu and Xu arguments has changed to Xr, Yr, Xu and Yu
- input type for the argument method has changed. Previously it received a 
character string  indicating the type of local regresion (i.e. "pls", 
"wapls1" or "gpr"). Now it receives an object of class local_fit which is output 
by the new local_fit fucntions. 
- dissimilarityM has been deprecated. It was used to pass a dissimilarity matrix 
computed outiside the mbl fucntion. This can be done now with the new argument 
diss_method of mbl which was previosly named "sm" and it was in mblControl()


### neigCleaning (now search_neighbors)
- Function renamed to search_neighbors
- Argument cores is deprecated, it wwas used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. Multi-threading via OpenMP is still 
available, but the number of threads must be set in the R global enviroment (e.g.
Sys.setenv("OMP_NUM_THREADS" = <a integer>)).
- Argument cores is deprecated, it wwas used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. Multi-threading via OpenMP is still 
available, but the number of threads must be set in the R global enviroment (e.g.
Sys.setenv("OMP_NUM_THREADS" = <a integer>)).


### mblControl (renamed to mbl_control):
- sm argument is deprecated. Now the dissmilarity metric is an argument of the 
mbl fucntion
- scale and center arguments have been moved to the mbl fucntion
- Argument range.pred.lim renamed to range_prediction_limits
- Argument cores is deprecated, it wwas used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. Multi-threading via OpenMP is still 
available, but the number of threads must be set in the R global enviroment (e.g.
Sys.setenv("OMP_NUM_THREADS" = <a integer>)).
- k0, pcMethod, ghMethod are deprecated
- localOptimization has been renamed to tune_locally
- valMethod has been renamed to validation_type
- Option "loc_crossval" in validation has been renamed to "local_cv"

### plot.mbl
- option "pca" was replaced by option "gh" which plots the pls projection used 
for computing the gh distance in mbl 


### Note

In `R`, the number of threads can be set by

Sys.setenv("OMP_NUM_THREADS" = 6)
