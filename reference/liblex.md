# Build a precomputed library of localised experts using memory-based learning

Constructs a library of local predictive models based on memory-based
learning (MBL). For each anchor observation, a local regression model is
fitted using its nearest neighbors from the reference set. This
implementation is based on the methods proposed in Ramirez-Lopez et al.
(2026b).

## Usage

``` r
liblex(Xr, Yr, neighbors,
       diss_method = diss_pca(ncomp = ncomp_by_opc()),
       fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
       anchor_indices = NULL, gh = TRUE, group = NULL,
       control = liblex_control(), verbose = TRUE, ...)

# S3 method for class 'liblex'
predict(object, newdata, diss_method = NULL,
        weighting = c("gaussian", "tricube", "triweight", "triangular",
                      "quartic", "parabolic", "cauchy", "none"),
        adaptive_bandwidth = TRUE, reliability_weighting = TRUE,
        range_prediction_limits = FALSE, residual_cutoff = NULL,
        enforce_indices = NULL, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        verbose = TRUE, allow_parallel = TRUE, blas_threads = 1L, ...)
        
# S3 method for class 'liblex'
plot(x, ...)
```

## Arguments

- Xr:

  A numeric matrix of predictor variables with dimensions `n × p`
  (observations in rows, variables in columns). Column names are
  required.

- Yr:

  A numeric vector or single-column matrix of length `n` containing the
  response variable values corresponding to `Xr`. Missing values are
  allowed (see Details).

- neighbors:

  A neighbor selection object specifying how to select neighbors. Use
  [`neighbors_k()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
  for fixed k-nearest neighbors or
  [`neighbors_diss()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
  for dissimilarity threshold-based selection.

- diss_method:

  For `liblex`: either a `diss_*` object specifying the dissimilarity
  method, or a precomputed numeric dissimilarity matrix. If a `diss_*`
  object (e.g.,
  [`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md),
  [`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md),
  [`diss_correlation()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md),
  [`diss_euclidean()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md)),
  dissimilarities are computed internally. Additional parameters
  (centering, scaling, number of components) are controlled within the
  object. If a matrix is provided:

  - When `anchor_indices = NULL`: must be square (`n × n`) with zeros on
    the diagonal.

  - When `anchor_indices` is specified (length `m`): must have
    dimensions `n × m`, with `diss_method[anchor_indices, ]` having
    zeros on the diagonal.

  Default is `diss_pca(ncomp = ncomp_by_opc())`.

  For `predict.liblex`: a dissimilarity method object created by
  `diss_*()` constructors. If not provided, uses the method stored in
  `object`. Required if `object` was built with a precomputed
  dissimilarity matrix.

- fit_method:

  A `local_fit` object specifying the local regression method. Currently
  supported:
  [`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  and
  [`fit_wapls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md).
  [`fit_gpr()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  is not yet supported. Default is
  `fit_wapls(min_ncomp = 3, max_ncomp = 15)`.

- anchor_indices:

  An optional integer vector specifying row indices of `Xr` around which
  to build local models. If `NULL` (default), models are built for all
  observations. See Details.

- gh:

  Logical indicating whether to compute the GH distance (Mahalanobis
  distance in PLS score space) for each anchor observation. Default is
  `TRUE`.

- group:

  An optional factor assigning group labels to observations in `Xr`.
  Used for leave-group-out validation to avoid pseudo-replication. When
  an observation is selected for validation, all observations from the
  same group are excluded from model fitting. Default is `NULL` (each
  observation is its own group).

- control:

  A list of control parameters created by
  [`liblex_control()`](https://l-ramirez-lopez.github.io/resemble/reference/liblex_control.md).
  Default is
  [`liblex_control()`](https://l-ramirez-lopez.github.io/resemble/reference/liblex_control.md).

- verbose:

  Logical indicating whether to display progress messages. Default is
  `TRUE`.

- ...:

  Additional arguments (currently unused).

- object:

  A fitted object of class `"liblex"`, created by `liblex` (for
  `predict`).

- newdata:

  A numeric matrix or data frame containing new predictor values. Must
  include all predictors used in `object`.

- weighting:

  Character string specifying the kernel weighting function applied to
  neighbours when combining predictions. Options are: `"gaussian"`
  (default), `"tricube"`, `"triweight"`, `"triangular"`, `"quartic"`,
  `"parabolic"`, `"cauchy"`, or `"none"` (equal weights). See Details
  for kernel definitions.

- adaptive_bandwidth:

  Logical indicating whether to use adaptive bandwidth for kernel
  weighting. When `TRUE` (default), the bandwidth is set to the maximum
  dissimilarity within the neighborhood of each observation, so weights
  adapt to local density. When `FALSE`, a fixed global bandwidth is used
  across all predictions, which may result in uniform weights in sparse
  regions or overly concentrated weights in dense regions.

- reliability_weighting:

  Logical indicating whether to weight expert predictions by their
  estimated reliability. When `TRUE` (default), the contribution of each
  experrt is additionally weighted by the inverse of its
  cross-validation residual variance, giving more influence to models
  that performed well during fitting. When `FALSE`, only
  dissimilarity-based kernel weights are used.

- probs:

  A numeric vector of probabilities in \\\[0, 1\]\\ for computing
  weighted quantiles of expert predictions. Default is
  `c(0.05, 0.25, 0.5, 0.75, 0.95)`.

- range_prediction_limits:

  Logical. If `TRUE`, predictions falling outside the 5th–95th
  percentile range of neighbour response values are clipped to those
  limits. Default is `FALSE`.

- residual_cutoff:

  Numeric threshold for excluding models. Models with absolute residuals
  exceeding this value are penalized during neighbour selection. Default
  is `NULL` (no exclusion).

- enforce_indices:

  Optional integer vector specifying model indices that must always be
  included in each prediction neighborhood. These models are assigned
  the minimum dissimilarity of the neighborhood to ensure selection.
  Default is `NULL` (no enforced models).

- allow_parallel:

  Logical indicating whether parallel computation is permitted if a
  backend is registered. Default is `TRUE`.

- blas_threads:

  Integer specifying the number of BLAS threads to use. Default is `1L`.
  Requires the RhpcBLASctl package for thread control.

- x:

  An object of class `"liblex"` as returned by `liblex()`.

## Value

**For `liblex`:** A list of class `"liblex"` (when
`control$mode = "build"`) or `"liblex_validation"` (when
`control$mode = "validate"`) containing:

- `dissimilarity`: List containing the dissimilarity method and matrix.

- `fit_method`: Fit constructor from `fit_method`.

- `gh`: If `gh = TRUE`, a list with GH distances and the PLS projection.

- `results`: Data frame of validation statistics for each parameter
  combination (if validation was performed).

- `best`: The optimal parameter combination based on `control$metric`.

- `optimal_params`: List with optimal `k` and `ncomp` values.

- `residuals`: Residuals from predictions using optimal parameters.

- `coefficients`: (Build mode only) List of regression coefficients:
  `B0` (intercepts), `B` (slopes).

- `vips`: (Build mode only) Variable importance in projection scores.

- `selectivity_ratios`: (Build mode only) Selectivity ratios for each
  predictor.

- `scaling`: (Build mode only) Centering and scaling vectors for
  prediction.

- `neighborhood_stats`: Statistics (response quantiles) for each
  neighborhood size.

- `anchor_indices`: The anchor indices used.

- `neighbors`: The object passed to `neighbors`.

**For `predict.liblex`:** A list with the following components:

- `predictions`: A data frame containing:

  - `pred`: Weighted mean predictions.

  - `pred_sd`: Weighted standard deviation of expert predictions.

  - `q*`: Weighted quantiles at probabilities specified by `probs`.

  - `gh`: Global Mahalanobis distance (if computed during fitting).

  - `min_yr`: Minimum response value (5th percentile) across neighbours.

  - `max_yr`: Maximum response value (95th percentile) across
    neighbours.

  - `below_min`: Logical indicating prediction below `min_yr`.

  - `above_max`: Logical indicating prediction above `max_yr`.

- `neighbors`: A list with:

  - `indices`: Matrix of neighbour indices (models) for each
    observation.

  - `dissimilarities`: Matrix of corresponding dissimilarity scores.

- `expert_predictions`: A list with:

  - `weights`: Matrix of kernel weights applied to each expert.

  - `predictions`: Matrix of raw predictions from each expert.

  - `weighted`: Matrix of weighted predictions from each expert.

## Details

By default, local models are constructed for all `n` observations in the
reference set. Alternatively, specify a subset of `m` observations
(`m < n`) via `anchor_indices` to reduce computation.

Each local model uses neighbors selected from the full reference set,
but models are only built for anchor observations. This is useful for
large datasets where building models for all observations is
computationally prohibitive.

When dissimilarity methods depend on `Yr` (e.g., PLS-based distances),
the response values of anchor observations are excluded during
dissimilarity computation for efficiency. However, anchor response
values are always used when fitting local models.

The number of anchors must not exceed 90% of `nrow(Xr)`; to build models
for all observations, use `anchor_indices = NULL`.

### Relationship between anchors and neighborhood size

The `neighbors` argument controls the neighborhood size (`k`) used both
for fitting local models and for retrieving experts during prediction.
When `anchor_indices` is specified, the number of available experts
equals the number of anchors. If `max(k)` exceeds the number of anchors
and tuning selects a large optimal `k`, prediction will retrieve fewer
experts than specified. For reliable predictions, ensure the number of
anchors is at least as large as the maximum `k` value being evaluated.

### Missing values in Yr

Missing values in `Yr` are permitted. Observations with missing response
values can still serve as neighbors but are excluded from model fitting
as target observations.

### GH distance

The GH distance is computed independently from `diss_method` using a PLS
projection with optimized component selection. This provides a measure
of how far each observation lies from the center of the reference set in
the PLS score space.

### Validation and tuning

When `control$mode = "validate"` or `control$tune = TRUE`,
nearest-neighbor cross-validation is performed. For each anchor
observation, its nearest neighbor is excluded, a model is fitted on
remaining neighbors, and the excluded neighbor's response is predicted.
This provides validation statistics for parameter selection.

### Prediction

For each observation in `newdata`, the `predict` method:

1.  Computes dissimilarities to anchor observations (or their
    neighbourhood centres) stored in `object`.

2.  Selects the `k` nearest neighbours based on the optimal `k`
    determined during model fitting.

3.  Applies kernel weighting based on dissimilarity.

4.  Combines expert predictions using weighted averaging.

### Kernel weighting functions

The weighting functions follow Cleveland and Devlin (1988). Let \\d\\ be
the normalised dissimilarity (scaled to \\\[0, 1\]\\ within the
neighbourhood when `adaptive_bandwidth = TRUE`). The available kernels
are:

- `"gaussian"`: \\w = \exp(-d^2)\\

- `"tricube"`: \\w = (1 - d^3)^3\\

- `"triweight"`: \\w = (1 - d^2)^3\\

- `"triangular"`: \\w = 1 - d\\

- `"quartic"`: \\w = (1 - d^2)^2\\

- `"parabolic"`: \\w = 1 - d^2\\

- `"cauchy"`: \\w = 1 / (1 + d^2)\\

- `"none"`: \\w = 1\\ (equal weights)

## References

Cleveland, W. S., & Devlin, S. J. (1988). Locally weighted regression:
An approach to regression analysis by local fitting. *Journal of the
American Statistical Association*, 83(403), 596–610.

Naes, T., Isaksson, T., & Kowalski, B. (1990). Locally weighted
regression and scatter correction for near-infrared reflectance data.
*Analytical Chemistry*, 62(7), 664–673.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. (2013). The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex datasets.
Geoderma, 195-196, 268-279.

Ramirez-Lopez, L., Metz, M., Lesnoff, M., Orellano, C., Perez-Fernandez,
E., Plans, M., Breure, T., Behrens, T., Viscarra Rossel, R., & Peng, Y.
(2026b). Rethinking local spectral modelling: From per-query refitting
to model libraries. *Analytica Chimica Acta*, under review.

Rajalahti, T., Arneberg, R., Berven, F.S., Myhr, K.M., Ulvik, R.J.,
Kvalheim, O.M. (2009). Biomarker discovery in mass spectral profiles by
means of selectivity ratio plot. Chemometrics and Intelligent Laboratory
Systems, 95(1), 35-48.

## See also

[`liblex_control()`](https://l-ramirez-lopez.github.io/resemble/reference/liblex_control.md)
for control parameters,
[`neighbors_k()`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)
for neighborhood specification,
[`diss_pca()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pca.md),
[`diss_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_pls.md),
[`diss_correlation()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_correlation.md)
for dissimilarity methods,
[`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md),
[`fit_wapls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
for fitting methods.

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
if (FALSE) { # \dontrun{
library(prospectr)
data(NIRsoil)

# Preprocess spectra
NIRsoil$spc_pr <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
  m = 1, p = 1, w = 7
)

# Missing values in the response are allowed
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]

# Build library
model_library <- liblex(
  Xr = train_x,
  Yr = train_y,
  neighbors = neighbors_k(c(30, 40)),
  diss_method = diss_correlation(ws = 27, scale = TRUE),
  fit_method = fit_wapls(
    min_ncomp = 4,
    max_ncomp = 17,
    scale = FALSE,
    method = "mpls"
  ),
  control = liblex_control(tune = TRUE)
)

# Visualise neighborhood centroids and samples to predict
matplot(
  as.numeric(colnames(model_library$scaling$local_x_center)),
  t(test_x),
  col = rgb(1, 0, 0, 0.3),
  lty = 1,
  type = "l",
  xlab = "Wavelength (nm)",
  ylab = "First derivative detrended absorbance"
)
matlines(
  as.numeric(colnames(model_library$scaling$local_x_center)),
  t(model_library$scaling$local_x_center),
  col = rgb(0, 0, 1, 0.3),
  lty = 1,
  type = "l"
)
grid(lty = 1)
legend(
  "topright",
  legend = c("Samples to predict", "Neighborhood centroids"),
  col = c(rgb(1, 0, 0, 0.8), rgb(0, 0, 1, 0.8)),
  lty = 1,
  lwd = 2,
  bty = "n"
)

# Predict new observations
y_hat_liblex <- predict(model_library, test_x)

# Predicted versus observed values
lims <- range(y_hat_liblex$predictions$pred, test_y, na.rm = TRUE)
plot(
  y_hat_liblex$predictions$pred,
  test_y,
  pch = 16,
  col = rgb(0, 0, 0, 0.5),
  xlab = "Predicted",
  ylab = "Observed",
  xlim = lims,
  ylim = lims
)
abline(a = 0, b = 1, col = "red")
grid(lty = 1)

## run liblex in parallel (requires a parallel backend, e.g., doParallel)
library(doParallel)
n_cores <- min(2, parallel::detectCores())
clust <- makeCluster(n_cores)
registerDoParallel(clust)

model_library2 <- liblex(
  Xr = train_x,
  Yr = train_y,
  neighbors = neighbors_k(c(30, 40)),
  fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, method = "simpls")
)

y_hat_liblex2 <- predict(model_library2, test_x)
registerDoSEQ()
try(stopCluster(clust))
} # }
```
