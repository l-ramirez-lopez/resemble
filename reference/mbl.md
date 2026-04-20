# Memory-based learning (mbl)

Memory-based learning (a.k.a. instance-based learning or local
regression) is a non-linear lazy learning approach for predicting a
response variable from predictor variables. For each observation in a
prediction set, a local regression is fitted using a subset of similar
observations (nearest neighbors) from a reference set. This function
does not produce a global model.

## Usage

``` r
mbl(Xr, Yr, Xu, Yu = NULL,
    neighbors,
    diss_method = diss_pca(ncomp = ncomp_by_opc()),
    diss_usage = c("none", "predictors", "weights"),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
    spike = NULL, group = NULL,
    gh = FALSE,
    control = mbl_control(),
    verbose = TRUE, seed = NULL,
    k, k_diss, k_range, method, pc_selection,
    center, scale, documentation, ...)

# S3 method for class 'mbl'
plot(x, what = c("validation", "gh"), metric = "rmse", ncomp = c(1, 2), ...)

get_predictions(x)

# S3 method for class 'mbl'
plot(x, what = c("validation", "gh"), metric = "rmse", ncomp = c(1, 2), ...)
```

## Arguments

- Xr:

  A matrix of predictor variables for the reference data (observations
  in rows, variables in columns). Column names are required.

- Yr:

  A numeric vector or single-column matrix of response values
  corresponding to `Xr`. NA values are not permitted.

- Xu:

  A matrix of predictor variables for the data to be predicted
  (observations in rows, variables in columns). Must have the same
  column names as `Xr`.

- Yu:

  An optional numeric vector or single-column matrix of response values
  corresponding to `Xu`. Used for computing prediction statistics.
  Default is `NULL`.

- neighbors:

  A neighbor selection object specifying how to select neighbors. Use
  [`neighbors_k()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
  for fixed k-nearest neighbors or
  [`neighbors_diss()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
  for dissimilarity threshold-based selection.

- diss_method:

  A dissimilarity method object or a precomputed dissimilarity matrix.
  Available constructors:

  - [`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md):
    Mahalanobis distance in PCA score space. This is the default where
    the number of components is optimized using side information (see
    `ncomp_by_opc`()).

  - [`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md):
    Mahalanobis distance in PLS score space

  - [`diss_euclidean()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md):
    Euclidean distance

  - [`diss_mahalanobis()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md):
    Mahalanobis distance

  - [`diss_cosine()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md):
    Cosine dissimilarity

  - [`diss_correlation()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md):
    Correlation-based dissimilarity

  A precomputed matrix can also be passed. When
  `diss_usage = "predictors"`, it must be square with dimensions
  `(nrow(Xr) + nrow(Xu))` and zeros on the diagonal. Otherwise, it must
  have `nrow(Xr)` rows and `nrow(Xu)` columns.

- diss_usage:

  How dissimilarity information is used in local models:

  - `"none"` (default): dissimilarities used only for neighbor selection

  - `"predictors"`: local dissimilarity matrix columns added as
    predictors

  - `"weights"`: neighbors weighted by dissimilarity using a tricubic
    function

- fit_method:

  A local fitting method object. Available constructors:

  - [`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md):
    Partial least squares regression

  - [`fit_wapls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md):
    Weighted average PLS (default)

  - [`fit_gpr()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md):
    Gaussian process regression

- spike:

  An integer vector indicating indices of observations in `Xr` to force
  into (positive values) or exclude from (negative values) all
  neighborhoods. Default is `NULL`. Spiking does not change neighborhood
  size; forced observations displace the most distant neighbors.

- group:

  An optional factor assigning group labels to `Xr` observations (e.g.,
  measurement batches). Used to avoid pseudo-replication in
  cross-validation: when one observation is held out, all observations
  from its group are also removed.

- gh:

  Logical indicating whether to compute global Mahalanobis (GH)
  distances. Default is `FALSE`. GH distances measure how far each
  observation lies from the center of the reference set in PLS score
  space. The computation uses a fixed methodology: PLS projection with
  the number of components selected via
  [`ncomp_by_opc()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)
  (capped at 40). This is independent of the `diss_method` argument.

- control:

  A list from
  [`mbl_control()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl_control.md)
  specifying validation type, tuning options, and other settings.

- verbose:

  Logical indicating whether to display a progress bar. Default is
  `TRUE`. Not shown during parallel execution.

- seed:

  An integer for random number generation, enabling reproducible
  cross-validation results. Default is `NULL`.

- k:

  Deprecated.

- k_diss:

  Deprecated.

- k_range:

  Deprecated.

- method:

  Deprecated.

- pc_selection:

  Deprecated.

- center:

  Deprecated.

- scale:

  Deprecated.

- documentation:

  Deprecated.

- ...:

  Additional arguments (currently unused).

- x:

  An object of class `mbl` (as returned by `mbl`).

- what:

  Character vector specifying what to plot. Options are `"validation"`
  (validation statistics) and/or `"gh"` (PLS scores used for GH distance
  computation). Default is both.

- metric:

  Character string specifying which validation statistic to plot.
  Options are `"rmse"`, `"st_rmse"`, or `"r2"`. Only used when
  `"validation"` is in `what`.

- ncomp:

  Integer vector of length 1 or 2 specifying which PLS components to
  plot. Default is `c(1, 2)`. Only used when `"gh"` is in `what`.

## Value

### mbl

For `mbl()`, a list of class `mbl` containing:

- `control`: control parameters from `control`

- `fit_method`: fit constructor from `fit_method`

- `Xu_neighbors`: list with neighbor indices and dissimilarities

- `dissimilarities`: dissimilarity method and matrix (if
  `return_dissimilarity = TRUE` in `control`)

- `n_predictions`: number of predictions made

- `gh`: GH distances for `Xr` and `Xu` (if `gh = TRUE`)

- `validation_results`: validation statistics by method

- `results`: list of data.frame objects with predictions, one per
  neighborhood size

- `seed`: the seed value used

Each results table contains:

- `o_index`: observation index

- `k`: number of neighbors used

- `k_diss`, `k_original`: (`neighbors_diss` only) threshold and original
  count

- `ncomp`: (`fit_pls` only) number of PLS components

- `min_ncomp`, `max_ncomp`: (`fit_wapls` only) component range

- `yu_obs`, `pred`: observed and predicted values

- `yr_min_obs`, `yr_max_obs`: response range in neighborhood

- `index_nearest_in_Xr`, `index_farthest_in_Xr`: neighbor indices

- `y_nearest`, `y_farthest`: neighbor response values

- `diss_nearest`, `diss_farthest`: neighbor dissimilarities

- `y_nearest_pred`: (NNv validation) leave-one-out prediction

- `loc_rmse_cv`, `loc_st_rmse_cv`: (local_cv validation) CV statistics

- `loc_ncomp`: (local dissimilarity only) components used locally

### Get predictions

The `get_predictions()` function extracts predicted values from an
object of class `mbl`. It returns a `data.frame` containing the
predictions.

## Details

### Spiking

The `spike` argument forces specific reference observations into or out
of neighborhoods. Positive indices are always included; negative indices
are always excluded. When observations are forced in, the most distant
neighbors are displaced to maintain neighborhood size. See Guerrero et
al. (2010).

### Dissimilarity usage

When `diss_usage = "predictors"`, the local dissimilarity matrix columns
are appended as additional predictor variables, which can improve
predictions (Ramirez-Lopez et al., 2013a).

When `diss_usage = "weights"`, neighbors are weighted using a tricubic
function (Cleveland and Devlin, 1988; Naes et al., 1990):

\\W\_{j} = (1 - v^{3})^{3}\\

where \\v = d(xr_i, xu_j) / \max(d)\\.

### GH distance

The global Mahalanobis distance (GH) measures how far each observation
lies from the center of the reference set. It is always computed using a
PLS projection with the number of components optimized via
[`ncomp_by_opc()`](https://l-ramirez-lopez.github.io/resemble/reference/ncomp_selection.md)
(maximum 40 components or `nrow(Xr)`, whichever is smaller). This
methodology is fixed and independent of the `diss_method` specified for
neighbor selection.

GH distances are useful for identifying extrapolation: observations with
high GH values lie far from the calibration space and may yield
unreliable predictions.

### Grouping

The `group` argument enables leave-group-out cross-validation. When
`validation_type = "local_cv"` in
[`mbl_control()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl_control.md),
the `p` parameter refers to the proportion of groups (not observations)
retained per iteration.

## References

Cleveland, W. S., and Devlin, S. J. 1988. Locally weighted regression:
an approach to regression analysis by local fitting. Journal of the
American Statistical Association 83:596-610.

Guerrero, C., Zornoza, R., Gomez, I., Mataix-Beneyto, J. 2010. Spiking
of NIR regional models using observations from target sites: Effect of
model size on prediction accuracy. Geoderma 158:66-77.

Naes, T., Isaksson, T., Kowalski, B. 1990. Locally weighted regression
and scatter correction for near-infrared reflectance data. Analytical
Chemistry 62:664-673.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex data sets.
Geoderma 195-196:268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R.,
Dematte, J.A.M., Scholten, T. 2013b. Distance and similarity-search
metrics for use with soil vis-NIR spectra. Geoderma 199:43-53.

Rasmussen, C.E., Williams, C.K. 2006. Gaussian Processes for Machine
Learning. MIT Press.

Shenk, J., Westerhaus, M., Berzaghi, P. 1997. Investigation of a LOCAL
calibration procedure for near infrared instruments. Journal of Near
Infrared Spectroscopy 5:223-232.

## See also

[`mbl_control`](https://l-ramirez-lopez.github.io/resemble/reference/mbl_control.md),
[`neighbors_k`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md),
[`neighbors_diss`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md),
[`diss_pca`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md),
[`diss_pls`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md),
[`fit_pls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md),
[`fit_wapls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md),
[`fit_gpr`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md),
[`search_neighbors`](https://l-ramirez-lopez.github.io/resemble/reference/search_neighbors.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120) and
Antoine Stevens

## Examples

``` r
if (FALSE) { # \dontrun{
library(prospectr)
data(NIRsoil)

# Preprocess: detrend + first derivative with Savitzky-Golay
sg_det <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
  m = 1, p = 1, w = 7
)
NIRsoil$spc_pr <- sg_det

# Split data
test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$CEC), ]
test_y <- NIRsoil$CEC[NIRsoil$train == 0 & !is.na(NIRsoil$CEC)]
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$CEC), ]
train_y <- NIRsoil$CEC[NIRsoil$train == 1 & !is.na(NIRsoil$CEC)]

# Example 1: Spectrum-based learner (Ramirez-Lopez et al., 2013)
ctrl <- mbl_control(validation_type = "NNv")

sbl <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  neighbors = neighbors_k(seq(40, 140, by = 20)),
  diss_method = diss_pca(ncomp = ncomp_by_opc(40)),
  fit_method = fit_gpr(),
  control = ctrl
)
sbl
plot(sbl)
get_predictions(sbl)

# Example 2: With known Yu
sbl_2 <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  neighbors = neighbors_k(seq(40, 140, by = 20)),
  fit_method = fit_gpr(),
  control = ctrl
)
plot(sbl_2)

# Example 3: LOCAL algorithm (Shenk et al., 1997)
local_algo <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  neighbors = neighbors_k(seq(40, 140, by = 20)),
  diss_method = diss_correlation(),
  diss_usage = "none",
  fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
  control = ctrl
)
plot(local_algo)

# Example 4: Using dissimilarity as predictors
local_algo_2 <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  neighbors = neighbors_k(seq(40, 140, by = 20)),
  diss_method = diss_pca(ncomp = ncomp_by_opc(40)),
  diss_usage = "predictors",
  fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
  control = ctrl
)
plot(local_algo_2)

# Example 5: Parallel execution
library(doParallel)
n_cores <- min(2, parallel::detectCores())
clust <- makeCluster(n_cores)
registerDoParallel(clust)

local_algo_par <- mbl(
  Xr = train_x,
  Yr = train_y,
  Xu = test_x,
  Yu = test_y,
  neighbors = neighbors_k(seq(40, 140, by = 20)),
  diss_method = diss_correlation(),
  fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
  control = ctrl
)

registerDoSEQ()
try(stopCluster(clust))
} # }
```
