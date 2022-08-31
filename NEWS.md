# `resemble`


`resemble 2.2.1 (Fix-Hodges)`
===============

### Improvements and fixes

* Fixed: An error was thrown when passing a pre-computed distance matrix to the 
`diss_method` argument in `mbl()` ([#24](https://github.com/l-ramirez-lopez/resemble/issues/24)).

* Documentation is now compatible with HTML5.


`resemble 2.1 (piapia)`
===============

23.09.2021

### New features

* The argument `seed` was added to the mbl function. It is used to gurantee 
reproducibility of cross-validation results. 

* A modified PLS method was implemented (see `local_fit_pls()` and `local_fit_wapls()`).
It uses correlation bteween response and predictors to derive the PLS weights. 

### Improvements and fixes

* A Bug in the computation of the explained variance of X for pls models was 
detected and fixed. The pls related functions were underestimating the amount 
of variance explained by each PLS component. The explained variance was being 
computed from the matrix of scores, in this version it is computed from the 
reconstructed spectra at each PLS iteration.

* Manual selection of components in `pc_projection()` and `pls_projection()` now 
accepts to select only 1 component (before the minimum was 2).

* `pls_projection()` now includes a new argument (`method`) to allow the user to select 
between the standard pls algorithm and a modified pls algorithm.

* `ortho_diss()`, `dissimilarity()` and `search_neighbors()` functions include a new 
dissimilarity method: `"mpls"` (modified pls).

* An internal function for stratified random sampling for cross-validtaion 
purposes has been improved for computational speed. 

* The package was stripping some symbols for Rcpp functions in Makevars in order 
to reduce the installation size of the package. Now these lines have been 
commented to comply with CRAN policies.


`resemble 2.0 (gordillo)`
===============

02.07.2020 

During the recent lockdown we had the chance to inevest a enough time on the 
development of a new version of the package `resemble`. This new version comes 
with significant improvements as well as new functionality. For example, it now 
matches the tidyverse R style guide, it includes unit tests, includes new 
functionality for modeling, mbl is faster and a less memory-hungry.

### New features

* `search_neighbors()` function.

* `dissimilarity()` function.

### Improvements and fixes

* `mbl` is faster and a less memory-hungry.

* New vignette.

* unit tests have been introduced using the testthat package.



### Breaking changes

#### `orthoProjection()`, `pcProjection()`, `plsProjection()` (renamed to `ortho_projection()`, 
`pc_projection()`, `pls_projection()` respectively):

* `X2` argument renamed to `Xu` (for consistency throughout all the functions).

* Argument `scaled` renamed to `.scale`.

* Argument `max.iter` renamed to `max_iter`.

* Bug fix: when the `"pca.nipals"`" method was used and the method to select the pcs was `"opc"`,
  the function was returning 1 component less than the maximum requested.
  
* `"pca.nipals"` is now implemented in C++ via Rcpp.

* Bug fix in `plsProjection()`: when `"cumvar"` was used as the `pcSelection` method, an 
  error about data allocation in a matrix was retrieved.
  
* Argument `pcSelection` renamed to `pc_selection`.

* `...` is deprecated in both `pcProjection()` and `plsProjection()`.

* Argument `cores` is deprecated, it was used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. 

* Names the following outputs have been changed: `X.loadings`, `Y.loadings`, `sc.sdv` 
and `n.components`. Their new names are: `X_loadings`, `Y_loadings`, `sc_sdv` and 
`n_components`.


#### `corDiss()` (renamed to `cor_diss()`):

* `X2` argument renamed to `Xu` (for consistency throughout all the functions).

* Argument `scaled` renamed to `.scale`.

* default for `.scale` has changed from `TRUE` to `FALSE`.

* the dimnames of the resulting matrix are now Xr_1... Xr_n (previusly Xr.1... Xr.n).


#### `fDiss()` (renamed to `f_diss()`):

* `X2` argument renamed to `Xu` (for consistency throughout all the functions).

* Argument `scaled` renamed to `.scale`.

* default for `.scale` has changed from `TRUE` to `FALSE`.

* the dimnames of the resulting matrix are now Xr_1... Xr_n (previusly Xr.1... Xr.n).

* argument method changed to diss_method.


#### `sid()`:

* `X2` argument renamed to `Xu` (for consistency throughout all the functions).

* Argument `scaled` renamed to `.scale`.

* default for `.scale` has changed from `TRUE` to `FALSE`.

* the dimnames of the resulting matrix are now Xr_1... Xr_n (previusly Xr.1... Xr.n).



#### orthoDiss (renamed to ortho_diss):

* `X2` argument renamed to `Xu` (for consistency throughout all the functions).

* Argument `scaled` renamed to `.scale`.

* Argument `local` renamed to `.local`.

* Argument `pcSelection` renamed to `pc_selection`.

* Argument `return.all` renamed to `compute_all`.

* Argument `cores` is deprecated, it wwas used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. 

* When `.local = TRUE` a new output is produced: `neighborhood_info` which 
is a data.frame containing the relevant information about the neighborhood of 
each sample (e.g. neighborhood indices, number of components used at each 
neighborhood, etc).

* Output `global.variance.info` has been renamed to `global_variance_info`


#### `simEval()` (renamed to `sim_eval()`):

* argument `sideInf` renamed to `side_info`.

* argument `lower.tri` renamed to `lower_triangle`.

* argument `cores` renamed to `omp_threads`.

* `lower_triangle` is deprecated. Now if a vector is passed to d, the function assumes
  that it is a lower triangle of a distance matrix.
  
* Argument `cores` is deprecated, it was used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp.
  
  
#### `mbl()`

* `pls.max.iter`, `pls.tol` and `noise.v` were moved to `mbl()` from `mbl_control()`.

* Argument scaled (from `mbl_control()`) renamed to .scale and moved to `mbl()`.

* new arguments: `gh` and `spike`.

* Argument `pcSelection` renamed to `pc_selection`.

* Argument `mblCtrl` renamed to `control`.

* Argument `dissUsage` renamed to `diss_usage`.

* order of the `Yr`, `Xr`, `Yu` and `Xu` arguments has changed to `Xr`, `Yr`, `Xu` and `Yu`.

* input type for the argument method has changed. Previously it received a 
character string  indicating the type of local regresion (i.e. "pls", 
"wapls1" or "gpr"). Now it receives an object of class `local_fit` which is output 
by the new `local_fit` functions. 

* `dissimilarityM` has been deprecated. It was used to pass a dissimilarity matrix 
computed outside the `mbl()` function. This can be done now with the new argument 
`diss_method` of `mbl` which was previously named `"sm"` and it was in `mblControl()`.


#### `neigCleaning()` (now `search_neighbors()`)

* Function renamed to `search_neighbors`.

* Argument `cores` is deprecated, it was used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. 


#### `mblControl()` (renamed to `mbl_control()`):

* `sm` argument is deprecated. Now the dissmilarity metric is an argument of the 
mbl function.

* `scale` and `center` arguments have been moved to the `mbl()` function.

* Argument `range.pred.lim` renamed to `range_prediction_limits`.

* Argument `cores` is deprecated, it was used to set the number of cores in some 
c++ internal functions via OpenMP in Rcpp. 

* `k0`, `pcMethod`, `ghMethod` are deprecated.

* `localOptimization` has been renamed to `tune_locally`.

* `valMethod` has been renamed to `validation_type`.

* Option `"loc_crossval"` in validation_type has been renamed to `"local_cv"`.

#### `plot.mbl()`

* option `"pca"` was replaced by option `"gh"` which plots the pls projection used 
for computing the gh distance in `mbl()`. 



`resemble 1.3` (never released)
===============

11.04.2020 

The option 'movcor' for the argument sm of mblControl() is deprecated. The 
'movcor' moving window correlations between spectra as dissimilarity measure. 
Now This option can be specified by using 'cor' as the method in the argument 
'sm' and passing a window size to the argument 'ws'of mblControl(). If 'ws' 
is not specified, the standard correlation between spectra is computed.  
             
27.02.2020 

The argument 'resampling' in mblControl() has been renamed to 'number'
             
18.07.2019 

A bug in the scaling of the euclidean distances in fDiss was detected and fixed. 
The distance ratios (between samples) were correctly calculated, but the final 
scaling of the results was not properly done. The distance between Xi and Xj 
were scaled by taking the squared root of the mean of the squared differences 
and dividing it by the number of variables i.e. sqrt(mean((Xi-Xj)^2))/ncol(Xi),
however the correct calculation is done by taking the mean of the squared 
differences, dividing it by the number of variables and then compute the squared 
root i.e. sqrt(mean((Xi-Xj)^2)/ncol(Xi)). This bug had no effect on the 
computations of the nearest neighbors. 

`resemble 1.2.0001 (alma de coco)`
===============

13.09.2016 

A bug in the computation of the Mahalanobis distance in the PLS space was fixed. 
             
06.09.2016 

Thanks to Matthieu Lesnoff who found a bug in the predict.orthoProjection 
function (an error was thrown when PCA preditions were requested). This bug has 
been fixed. 
             
10.08.2016 

A bug in plot.mbl was fixed. It was not possible to plot mbl results when the 
k.diss argument (threshold distances) was used in the mbl function.   

10.08.2016 

Since the previous release, the "wapls1" regression (in the mbl function) 
is actually compatible with valMethod = "loc_crossval" (in the mblControl). 
In the previous documentation was wrongly stated otherwise. Now this has been 
corrected in the documentation.

09.08.2016 

the projection Matrix (projectionM) returned by plsProjection now only contains 
the columns corresponding only to the number components retrieved
             
`resemble 1.2 (alma de coco)`
===============

04.03.2016 

A patch was released for and extrange bug that prevented to run mbl 
in parallel when the gpr method was used.

19.01.2016 

Now it is possible to locally optimize the maximum and minimum pls factors in 
wapls1 local regressions.

09.12.2015 

Many thanks to Eva Ampe and Lorenzo Menichetti for their suggestions.

08.12.2015 

A method for better estimates of RMSE values computed for the 'wapls1' method 
has been implemented.

08.12.2015 

The 'wapls2' method of the mbl function is no longer supported because of several 
drawbacks computing reliable uncertainty estimates.

18.11.2015 

Several functions are now based on C++ for faster computations.

23.04.2014 

Added default variable names when they are missing and an error message when the 
names of Xr do not match the names of Xu.

23.04.2014 

plot.mbl draws now the circles around the actual centre function when the 
spectra is not centred for mbl.

20.03.2013 

The function movcorDist was removed since. it was included by mistake in the 
previous version of the package. The corDiss function  can be used in 
raplacement of movcorDist.

`resemble 1.1.1`
===============

19.03.2013 Hello world! Initial release of the package
