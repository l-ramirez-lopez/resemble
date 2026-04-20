# Evolutionary sample search for context-specific calibrations

Implements an evolutionary search algorithm that selects a subset from
large reference datasets (e.g., spectral libraries) to build
context-specific calibrations. The algorithm iteratively removes weak or
non-informative samples based on prediction error, spectral
reconstruction error, or dissimilarity criteria. This implementation is
based on the methods proposed in Ramirez-Lopez et al. (2026a).

## Usage

``` r
# Default S3 method
gesearch(Xr, Yr, Xu, Yu = NULL, Yu_lims = NULL,
         k, b, retain = 0.95, target_size = k,
         fit_method = fit_pls(ncomp = 10),
         optimization = "reconstruction",
         group = NULL, control = gesearch_control(),
         intermediate_models = FALSE,
         verbose = TRUE, seed = NULL, pchunks = 1L, ...)

# S3 method for class 'formula'
gesearch(formula, train, test, k, b, target_size, fit_method,
         ..., na_action = na.pass)

# S3 method for class 'gesearch'
predict(object, newdata, type = "response",
         what = c("final", "all_generations"), ...)

# S3 method for class 'gesearch'
plot(x, which = c("weakness", "removed"), ...)
```

## Arguments

- Xr:

  A numeric matrix of predictor variables for the reference data
  (observations in rows, variables in columns).

- Yr:

  A numeric vector or single-column matrix of response values
  corresponding to `Xr`. Only one response variable is supported.

- Xu:

  A numeric matrix of predictor variables for target observations (same
  structure as `Xr`).

- Yu:

  An optional numeric vector or single-column matrix of response values
  for `Xu`. Required when `optimization` includes `"response"`. Default
  is `NULL`.

- Yu_lims:

  A numeric vector of length 2 specifying expected response limits for
  the target population. Used with `optimization = "range"`.

- k:

  An integer specifying the number of samples in each resampling subset
  (gene size).

- b:

  An integer specifying the target average number of times each training
  sample is evaluated per iteration. Higher values (e.g., \>40) produce
  more stable results but increase computation time.

- retain:

  A numeric value in (0, 1\] specifying the proportion of samples
  retained per iteration. Default is 0.95. Values \>0.9 are recommended
  for stability. See
  [`gesearch_control`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch_control.md)
  for retention strategy.

- target_size:

  An integer specifying the target number of selected samples (gene pool
  size). Must be \>= `k`. Default is `k`.

- fit_method:

  A fit method object created with
  [`fit_pls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md).
  Specifies the regression model and scaling used during the search.
  Currently only
  [`fit_pls()`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  is supported.

- optimization:

  A character vector specifying optimization criteria:

  - `"reconstruction"`: (default) Retains samples based on spectral
    reconstruction error of `Xu` in PLS space.

  - `"response"`: Retains samples based on RMSE of predicting `Yu`.
    Requires `Yu`.

  - `"similarity"`: Retains samples based on Mahalanobis distance
    between `Xu` and training samples in PLS score space.

  - `"range"`: Removes samples producing predictions outside `Yu_lims`.

  Multiple criteria can be combined, e.g.,
  `c("reconstruction", "similarity")`.

- group:

  An optional factor assigning group labels to training observations.
  Used for leave-group-out cross-validation to avoid pseudo-replication.

- control:

  A list created with
  [`gesearch_control`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch_control.md)
  containing additional algorithm parameters.

- intermediate_models:

  A logical indicating whether to store models for each intermediate
  generation. Default is `FALSE`.

- verbose:

  A logical indicating whether to print progress information. Default is
  `TRUE`.

- seed:

  An integer for random number generation to ensure reproducibility.
  Default is `NULL`.

- pchunks:

  An integer specifying the chunk size used for memory-efficient
  parallel processing. Larger values divide the workload into smaller
  pieces, which can help reduce memory pressure. Default is 1L.

- formula:

  A [`formula`](https://rdrr.io/r/stats/formula.html) defining the
  model.

- train:

  A data.frame containing training data with model variables.

- test:

  A data.frame containing test data with model variables.

- na_action:

  A function for handling missing values in training data. Default is
  [`na.pass`](https://rdrr.io/r/stats/na.fail.html).

- object:

  A fitted `gesearch` object (for `predict`).

- newdata:

  A matrix or data.frame of new observations. For formula-fitted models,
  a data.frame containing all predictor variables is accepted. For
  non-formula models, a matrix is required.

- type:

  A character string specifying the prediction type. Currently only
  `"response"` is supported.

- what:

  A character string specifying which models to use for prediction:
  `"final"` (default) for predictions from final models only, or
  `"all_generations"` for predictions from all intermediate generations
  plus the final models.

- x:

  A `gesearch` object (for `plot`).

- which:

  Character string specifying what to plot: `"weakness"` (maximum
  weakness scores per generation) or `"removed"` (cumulative samples
  removed).

- ...:

  Additional arguments passed to methods.

## Value

**For `gesearch`:** A list of class `"gesearch"` containing:

- `x_local`: Matrix of predictors for selected samples.

- `y_local`: Vector of responses for selected samples.

- `indices`: Indices of selected samples from original training set.

- `complete_iter`: Number of completed iterations.

- `iter_weakness`: List with iteration-level weakness statistics.

- `samples`: List of sample indices retained at each iteration.

- `n_removed`: data.frame of samples removed per iteration.

- `control`: Copy of control parameters.

- `fit_method`: Fit constructor from `fit_method`.

- `validation_results`: Cross-validation in the training only set
  validation on the test set using models built only with the samples
  found.

- `final_models`: Final PLS model containing coefficients, loadings,
  scores, VIP, and selectivity ratios.

- `intermediate_models`: List of models per generation (if
  `intermediate_models = TRUE`).

- `seed`: RNG seed used.

**For `predict.gesearch`:**

- If `what = "final"`: a prediction matrix with `nrow(newdata)` rows and
  one column per PLS component.

- If `what = "all_generations"`: a named list of generations, where each
  generation contains a prediction matrix as above.

## Details

The `gesearch` algorithm requires a large reference dataset (`Xr`) where
the sample search is conducted, target observations (`Xu`), and three
tuning parameters: `k`, `b`, and `retain`.

The target observations (`Xu`) should represent the population of
interest. These may be selected via algorithms like Kennard-Stone when
response values are unavailable.

The algorithm iteratively removes weak samples from `Xr` based on:

- Increased RMSE when predicting `Yu`

- Increased PLS reconstruction error on `Xu`

- Increased dissimilarity to `Xu` in PLS space

A resampling scheme identifies samples that consistently appear in
high-error subsets. These are labeled weak and removed. The process
continues until approximately `target_size` samples remain.

The `gesearch()` function also returns a final model fitted on the
selected samples, which can be used for prediction. This model is
internally validated by cross-validation using only the selected samples
from the training/reference set. If `Yu` is available, a model fitted
only on the selected reference samples is first used to predict the
target samples. The final model is then refitted using both the selected
reference samples and the target samples used to guide the search,
provided that response values are available for those target samples.

### Parameter guidance

- `k`: Number of samples per resampling subset. See Lobsey et al. (2017)
  for guidance.

- `b`: Resampling intensity. Higher values increase stability but
  computational cost.

- `retain`: Proportion retained per iteration. Values \>0.9 recommended.

### Prediction

The `predict` method generates predictions from a fitted `gesearch`
object. If the model was fitted with a formula, `newdata` is validated
and transformed to the appropriate model matrix.

When `what = "all_generations"`, the return value is a named list with
one element per generation, where each element contains a prediction
matrix. This option requires `intermediate_models = TRUE` during
fitting.

## References

Lobsey, C.R., Viscarra Rossel, R.A., Roudier, P., Hedley, C.B. 2017.
rs-local data-mines information from spectral libraries to improve local
calibrations. European Journal of Soil Science 68:840-852.

Kennard, R.W., Stone, L.A. 1969. Computer aided design of experiments.
Technometrics 11:137-148.

Rajalahti, T., Arneberg, R., Berven, F.S., Myhr, K.M., Ulvik, R.J.,
Kvalheim, O.M. 2009. Biomarker discovery in mass spectral profiles by
means of selectivity ratio plot. Chemometrics and Intelligent Laboratory
Systems 95:35-48.

Ramirez-Lopez, L., Viscarra Rossel, R., Behrens, T., Orellano, C.,
Perez-Fernandez, E., Kooijman, L., Wadoux, A. M. J.-C., Breure, T.,
Summerauer, L., Safanelli, J. L., & Plans, M. (2026a). When spectral
libraries are too complex to search: Evolutionary subset selection for
domain-adaptive calibration. *Analytica Chimica Acta*, under review.

## See also

[`fit_pls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md),
[`gesearch_control`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch_control.md),
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120), Claudio
Orellano, [Craig Lobsey](https://orcid.org/0000-0001-5416-4520),
[Raphael Viscarra Rossel](https://orcid.org/0000-0003-1540-4748)

## Examples

``` r
# \donttest{
library(prospectr)
data(NIRsoil)

# Preprocess
sg_det <- savitzkyGolay(
  detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
  m = 1, p = 1, w = 7
)
NIRsoil$spc_pr <- sg_det

# Split data
train_x <- NIRsoil$spc_pr[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso), ]
train_y <- NIRsoil$Ciso[NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)]
test_x <- NIRsoil$spc_pr[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso), ]
test_y <- NIRsoil$Ciso[NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)]

# Basic search with reconstruction and similarity optimizations
gs <- gesearch(
  Xr = train_x, Yr = train_y,
  Xu = test_x, Yu = test_y,
  k = 50, b = 100, retain = 0.97,
  target_size = 200,
  fit_method = fit_pls(ncomp = 15, method = "mpls"),
  optimization = c("reconstruction", "similarity"),
  control = gesearch_control(retain_by = "probability"),
  seed = 42
)
#> Generation 1: Initial genes (samples): 548 >>  | Individuals: 1096 >> Generation 2: Active genes (samples): 517 >>> | Individuals: 1034 >>> Generation 3: Active genes (samples): 485     | Individuals: 970     Generation 4: Active genes (samples): 457 >   | Individuals: 914 >  Generation 5: Active genes (samples): 429 >>  | Individuals: 858 >> Generation 6: Active genes (samples): 403 >>> | Individuals: 806 >>>Generation 7: Active genes (samples): 377     | Individuals: 754    Generation 8: Active genes (samples): 354 >   | Individuals: 708 >  Generation 9: Active genes (samples): 332 >>  | Individuals: 664 >> Generation 10: Active genes (samples): 312 >>> | Individuals: 624 >>>Generation 11: Active genes (samples): 292     | Individuals: 584    Generation 12: Active genes (samples): 274 >   | Individuals: 548 >  Generation 13: Active genes (samples): 256 >>  | Individuals: 512 >> Generation 14: Active genes (samples): 240 >>> | Individuals: 480 >>>Generation 15: Active genes (samples): 224     | Individuals: 448    Generation 16: Active genes (samples): 210 >   | Individuals: 420 >  
#> Fitting final model on 210 selected genes...
#> and 184 target genes with available response values...

# Predict
preds <- predict(gs, test_x)

# Plot progress
plot(gs)

plot(gs, which = "removed")


# With reconstruction and response optimization (requires Yu)
gs_response <- gesearch(
  Xr = train_x, Yr = train_y,
  Xu = test_x, Yu = test_y,
  k = 50, b = 100, retain = 0.97,
  target_size = 200,
  fit_method = fit_pls(ncomp = 15),
  optimization = c("reconstruction", "response"),
  seed = 42
)
#> Generation 1: Initial genes (samples): 548 >>  | Individuals: 1096 >> Generation 2: Active genes (samples): 517 >>> | Individuals: 1034 >>> Generation 3: Active genes (samples): 487     | Individuals: 974     Generation 4: Active genes (samples): 458 >   | Individuals: 916 >  Generation 5: Active genes (samples): 431 >>  | Individuals: 862 >> Generation 6: Active genes (samples): 406 >>> | Individuals: 812 >>>Generation 7: Active genes (samples): 384     | Individuals: 768    Generation 8: Active genes (samples): 361 >   | Individuals: 722 >  Generation 9: Active genes (samples): 339 >>  | Individuals: 678 >> Generation 10: Active genes (samples): 319 >>> | Individuals: 638 >>>Generation 11: Active genes (samples): 300     | Individuals: 600    Generation 12: Active genes (samples): 282 >   | Individuals: 564 >  Generation 13: Active genes (samples): 265 >>  | Individuals: 530 >> Generation 14: Active genes (samples): 250 >>> | Individuals: 500 >>>Generation 15: Active genes (samples): 235     | Individuals: 470    Generation 16: Active genes (samples): 219 >   | Individuals: 438 >  Generation 17: Active genes (samples): 208 >>  | Individuals: 416 >> 
#> Fitting final model on 208 selected genes...
#> and 184 target genes with available response values...

# Parallel processing
library(doParallel)
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
n_cores <- min(4, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

gs_parallel <- gesearch(
  Xr = train_x, Yr = train_y,
  Xu = test_x,
  k = 50, b = 100, retain = 0.97,
  target_size = 200,
  fit_method = fit_pls(ncomp = 15),
  pchunks = 3,
  seed = 42
)
#> Generation 1: Initial genes (samples): 548 >>  | Individuals: 1096 >> Generation 2: Active genes (samples): 531 >>> | Individuals: 1062 >>> Generation 3: Active genes (samples): 515     | Individuals: 1030    Generation 4: Active genes (samples): 499 >   | Individuals: 998 >   Generation 5: Active genes (samples): 484 >>  | Individuals: 968 >> Generation 6: Active genes (samples): 469 >>> | Individuals: 938 >>>Generation 7: Active genes (samples): 454     | Individuals: 908    Generation 8: Active genes (samples): 440 >   | Individuals: 880 >  Generation 9: Active genes (samples): 426 >>  | Individuals: 852 >> Generation 10: Active genes (samples): 413 >>> | Individuals: 826 >>>Generation 11: Active genes (samples): 400     | Individuals: 800    Generation 12: Active genes (samples): 388 >   | Individuals: 776 >  Generation 13: Active genes (samples): 376 >>  | Individuals: 752 >> Generation 14: Active genes (samples): 364 >>> | Individuals: 728 >>>Generation 15: Active genes (samples): 353     | Individuals: 706    Generation 16: Active genes (samples): 342 >   | Individuals: 684 >  Generation 17: Active genes (samples): 331 >>  | Individuals: 662 >> Generation 18: Active genes (samples): 321 >>> | Individuals: 642 >>>Generation 19: Active genes (samples): 311     | Individuals: 622    Generation 20: Active genes (samples): 301 >   | Individuals: 602 >  Generation 21: Active genes (samples): 291 >>  | Individuals: 582 >> Generation 22: Active genes (samples): 282 >>> | Individuals: 564 >>>Generation 23: Active genes (samples): 273     | Individuals: 546    Generation 24: Active genes (samples): 264 >   | Individuals: 528 >  Generation 25: Active genes (samples): 256 >>  | Individuals: 512 >> Generation 26: Active genes (samples): 248 >>> | Individuals: 496 >>>Generation 27: Active genes (samples): 240     | Individuals: 480    Generation 28: Active genes (samples): 232 >   | Individuals: 464 >  Generation 29: Active genes (samples): 225 >>  | Individuals: 450 >> Generation 30: Active genes (samples): 218 >>> | Individuals: 436 >>>Generation 31: Active genes (samples): 211     | Individuals: 422    Generation 32: Active genes (samples): 204 >   | Individuals: 408 >  
#> Fitting final model on 204 selected genes...

stopCluster(cl)
registerDoSEQ()
# }
```
