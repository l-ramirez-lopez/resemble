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
# \donttest{
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
#> Predicting...
#> 1/113__________                          2/113__________                          3/113__________                          4/113__________                          5/113__________                          6/113__________                          7/113__________                          8/113__________                          9/113__________                          10/113__________                          11/113__________                          12/113__________                          13/113__________                          14/113__________                          15/113__________                          16/113__________                          17/113__________                          18/113__________                          19/113__________                          20/113__________                          21/113__________                          22/113__________                          23/113__________                          24/113__________                          25/113__________                          26/113__________                          27/113__________                          28/113__________                          29/113__________                          30/113__________                          31/113__________                          32/113__________                          33/113__________                          34/113__________                          35/113__________                          36/113__________                          37/113__________                          38/113__________                          39/113__________                          40/113__________                          41/113__________                          42/113__________                          43/113__________                          44/113__________                          45/113__________                          46/113__________                          47/113__________                          48/113__________                          49/113__________                          50/113__________                          51/113__________                          52/113__________                          53/113__________                          54/113__________                          55/113__________                          56/113__________                          57/113__________                          58/113__________                          59/113__________                          60/113__________                          61/113__________                          62/113__________                          63/113__________                          64/113__________                          65/113__________                          66/113__________                          67/113__________                          68/113__________                          69/113__________                          70/113__________                          71/113__________                          72/113__________                          73/113__________                          74/113__________                          75/113__________                          76/113__________                          77/113__________                          78/113__________                          79/113__________                          80/113__________                          81/113__________                          82/113__________                          83/113__________                          84/113__________                          85/113__________                          86/113__________                          87/113__________                          88/113__________                          89/113__________                          90/113__________                          91/113__________                          92/113__________                          93/113__________                          94/113__________                          95/113__________                          96/113__________                          97/113__________                          98/113__________                          99/113__________                          100/113__________                          101/113__________                          102/113__________                          103/113__________                          104/113__________                          105/113__________                          106/113__________                          107/113__________                          108/113__________                          109/113__________                          110/113__________                          111/113__________                          112/113__________                          113/113__________
sbl
#> --- mbl predictions --- 
#> Predictions: 113 
#> _______________________________________________________ 
#> Dissimilarity 
#> Dissimilarity: PCA
#>   method            : pca 
#>   ncomp             :  
#>   center            : TRUE 
#>   scale             : FALSE 
#>   return_projection : FALSE 
#> _______________________________________________________ 
#> Local fit method 
#> Fitting method: gpr
#>    noise_variance : 0.001 
#>    center         : TRUE 
#>    scale          : TRUE 
#> _______________________________________________________ 
#> 
#> Nearest-neighbor validation 
#> 
#>    k rmse st_rmse    r2
#>   40 4.04  0.0924 0.543
#>   60 3.94  0.0903 0.527
#>   80 3.87  0.0886 0.561
#>  100 4.03  0.0921 0.547
#>  120 4.41  0.1010 0.525
#>  140 4.17  0.0953 0.575
#> _______________________________________________________ 
plot(sbl)
#> GH distance not available in this object.

get_predictions(sbl)
#>          k_40      k_60      k_80     k_100     k_120     k_140
#> 1    5.366291  7.857863  6.471958 11.553644 10.022977 10.807479
#> 2    1.400000  1.400000  2.244692  1.400000  1.400000  1.400000
#> 3    6.200000  6.200000  6.200000  6.200000  5.800000  5.800000
#> 4    1.500000  1.500000  1.500000  1.500000  1.000000  1.000000
#> 5   12.786351 16.639999 16.639999 16.639999  8.445008  3.400000
#> 6   10.741227  9.762140  9.945469  9.400751 10.149403 11.161453
#> 7    7.175268  7.934039  9.845900  8.875633  7.575546  8.636831
#> 8    9.491909 11.197229 10.067193 11.144713  9.256177 11.434195
#> 9    6.484870  6.200000  6.200000  6.200000  6.200000  9.794681
#> 10   3.990645  7.422981  5.442374  6.879462  7.723470  5.770113
#> 11   9.713616  9.252345  8.581060  8.038268  9.631575  9.621054
#> 12  10.435846  7.132317 11.795727 10.409522 11.289904  4.510632
#> 13  12.445042 12.653790 11.339653  5.912106  6.553059  5.564825
#> 14   9.008466  9.520396  9.299262  8.394983  8.976025  8.215273
#> 15   6.293976  5.800000  7.487466  4.731545  4.300000  5.920069
#> 16  13.325691 13.071587 14.201130 12.607336 14.144815 14.667272
#> 17   9.261265  8.356862  9.892288  7.843793  7.528563  8.454964
#> 18  10.662036 12.842824 14.173631 13.171822 12.929165 15.779389
#> 19   8.216828  8.037974  8.044760  8.549928  9.480975  8.368563
#> 20  12.104446 12.453960 12.183754 13.721782 14.358436 16.538139
#> 21   6.143411  6.909202  7.226713  7.039899  8.869267  8.625495
#> 22  10.000000 10.003606 10.415908 10.507091  9.884616 10.857010
#> 23  10.138272  9.517861  8.255953  7.200000  6.670000  5.800000
#> 24   8.521843  8.125711  8.583765  8.272176 10.579878 10.404383
#> 25  10.508257 10.173600 10.332051 10.736832 10.101550  7.886847
#> 26  11.638130 13.260668 14.030277 13.942299 12.951840 11.641618
#> 27  10.324976 10.310218  9.016014  7.502927  9.674001  8.756714
#> 28   7.868676  9.316037  9.631500  8.654352  8.808088  9.252244
#> 29   9.853316  9.713163  9.052354  9.048732  9.188598 10.559924
#> 30  11.653871  8.713296 14.885913 16.509260 16.416419 14.233668
#> 31  13.921732 14.915639 12.271151 13.379921 14.151391 14.316209
#> 32  13.275878 14.590961 13.979778 15.679856 14.682451 13.646024
#> 33  12.717554 12.187847 12.137354 11.698036 12.121562 11.159062
#> 34  10.378420  9.755125  9.233504  9.851568 10.073800  9.527818
#> 35   8.997707 10.566882 12.574610 14.512233 12.802054 15.958623
#> 36  13.574038 12.034658 12.355753 10.743981  9.968448  9.256924
#> 37   7.907652  9.182024  9.022714  9.553548  9.424902  8.767498
#> 38  17.733282 17.298007 16.891180 15.502363 13.603654 13.086604
#> 39  13.061451 15.282673 15.379359 12.795550 12.424361 14.411097
#> 40  15.022399 13.008745  5.674590  9.626950  6.113952 10.264368
#> 41  14.810784 15.430000 13.374140 14.315422 15.716648 13.496353
#> 42  15.349556 14.504917 15.009797 14.170153 14.107184 14.539702
#> 43  11.662214 13.194270 12.614112 12.955567 12.875969 12.296964
#> 44  10.107842  8.873966  7.914077  7.802091  9.252788  7.160030
#> 45   8.503117  7.727977  9.662640  7.599134 11.111832 10.566321
#> 46  12.600367 12.182612 12.213114 12.631466 13.382120 13.286391
#> 47   8.957335  9.491847  9.391779 11.090071  9.374947  7.986726
#> 48  12.977092 11.699254 10.516655  8.837146  7.154024  9.010066
#> 49   6.869890  3.715063  4.273051  9.753037 11.254509 13.343952
#> 50  11.313785 10.930896  9.823930 11.169202 11.537209 13.383434
#> 51  11.392450 11.596231 11.342576 10.952034 11.695078 11.228242
#> 52  17.860411 16.890976 23.293491 18.820457 15.772065 12.240745
#> 53  12.177207 12.224919 12.277866 12.772630 12.153088 12.676602
#> 54  11.500802 10.621304 10.579637 11.650999 11.748315  9.985794
#> 55  10.900330 12.822712 12.563540 11.462962 11.863601 11.207547
#> 56   7.766412  7.492879  7.676950  6.500000  4.343222  4.037606
#> 57   6.200000  6.200000  4.300000  1.000000  1.000000  1.000000
#> 58  12.698892 13.419558 16.176491 14.824755 15.284647 15.224164
#> 59  10.816789 11.548145 11.104443 10.921899 14.684989 17.976810
#> 60  13.361053 13.339632 14.615455 14.745720 13.985392 14.050010
#> 61  12.886497 13.261124 13.616530 13.808178 14.290226 15.052559
#> 62  11.162853 11.661043 11.143414  9.928962  9.842503 10.574709
#> 63  11.492103 10.518247 11.377401 11.561751 10.903662 12.585243
#> 64  10.392321 10.705730 10.160228  9.071279  7.409505 10.172324
#> 65  11.439273 10.125311 10.751921 11.374736 10.405722 10.584475
#> 66  14.900000 15.301247 14.082922 14.385987 14.264618 15.167233
#> 67  12.341220 13.116782 14.043505 14.877981 12.443572 14.462110
#> 68  16.078046 14.317179  7.051860  7.751403  8.606793 10.879094
#> 69  13.332371 13.237294 13.011266 12.403483 12.881230 13.137435
#> 70  14.192167 14.636157 14.032613 14.016655 14.321142 14.203532
#> 71  12.258692 11.486424 11.952244 13.160166 12.670883 12.485906
#> 72  13.346033 13.331636 12.908847 13.677977 13.951663 13.826711
#> 73  12.771731 12.926312 13.075881 13.229911 12.601887 12.202344
#> 74  16.250000 15.449386 13.698776 12.899793 13.052920 14.009891
#> 75  11.397405 12.641180 11.824772  9.437022 10.384917 10.740320
#> 76  13.226112 12.843773 14.389620 13.979302 14.118970 13.102062
#> 77  13.647171 14.177603 11.957214 12.100640 13.382339  9.688852
#> 78  14.436876 13.988524 14.877786 13.877694 13.239839 13.478638
#> 79  13.230695 12.768605 16.123626 13.873099 13.520403 17.589699
#> 80  16.418913 18.717802 21.608788 22.645064 22.140674 22.603690
#> 81  12.162707 12.296867 12.189581 10.587045  9.812929 12.718767
#> 82  15.366885 14.808405 13.871114 14.507603 15.756809 16.213608
#> 83  14.496038 14.899361 14.948215 13.684437 15.425627 15.439781
#> 84  16.250000 15.165450 16.250000 14.803230 15.401394 15.839782
#> 85  16.005582 15.926431 16.250000 15.657242 15.552129 14.565730
#> 86  14.115869 13.879103 14.191283 13.948457 13.926694 13.669235
#> 87  13.244613 11.549591 12.698795 11.425800  9.118248  8.013099
#> 88  16.226723 15.030083 13.319724 12.634075 13.145142 11.973696
#> 89  15.914164 16.250000 16.199239 15.706410 15.642499 14.711290
#> 90   9.403191 11.902494 13.691113 14.517586 12.315372 13.042562
#> 91  12.422347 14.625355 15.262279 14.969622 16.023407 16.500206
#> 92  16.377287 15.489727 15.388546 14.598580 15.192217 17.523964
#> 93  16.250000 16.108808 16.442944 16.746959 18.110373 20.629502
#> 94  15.067896 14.941900 14.463845 13.886122 12.018772 13.908500
#> 95  14.433821 14.782151 14.421120 14.339668 14.851204 13.493163
#> 96  18.184475 16.557100 16.340874 16.773579 17.293753 18.750429
#> 97  11.045572 11.371914 16.213532 16.455711 15.515433 13.625323
#> 98  15.065760 14.523183 13.539033 13.706641 15.477081 16.594760
#> 99  16.719329 16.436073 16.711166 18.324021 17.447449 18.376649
#> 100 14.946303 14.909773 13.367309 13.684833 13.091084 11.361289
#> 101 16.250000 15.652974 16.250000 17.414049 17.860136 16.848060
#> 102 16.038399 15.991752 17.137638 17.782422 18.296142 17.864616
#> 103 22.480624 22.374812 23.139418 22.719773 23.734362 20.890542
#> 104 15.497179 14.177751 13.343324 12.945959 12.549507 12.693844
#> 105 17.074721 17.288724 17.332191 16.776095 20.847244 22.012033
#> 106 18.581822 21.000000 17.383438 21.000000 19.846562 13.224230
#> 107 20.550664 22.844029 20.817240 24.236577 23.183144 24.323990
#> 108 20.600000 16.296913 21.205608 21.333459 17.446449 17.840069
#> 109 19.173913 20.508120 15.128486 22.531503 30.452686 31.402592
#> 110 33.846442 31.385648 32.243710 35.314639 36.886568 32.033060
#> 111 22.020193 25.590408 29.721540 31.417951 33.109691 29.929283
#> 112 19.666719 18.568568 19.514221 21.158401 21.986762 21.184338
#> 113 40.419765 41.731254 42.054347 41.337163 40.802153 42.984958

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
#> Predicting...
#> 1/113__________                          2/113__________                          3/113__________                          4/113__________                          5/113__________                          6/113__________                          7/113__________                          8/113__________                          9/113__________                          10/113__________                          11/113__________                          12/113__________                          13/113__________                          14/113__________                          15/113__________                          16/113__________                          17/113__________                          18/113__________                          19/113__________                          20/113__________                          21/113__________                          22/113__________                          23/113__________                          24/113__________                          25/113__________                          26/113__________                          27/113__________                          28/113__________                          29/113__________                          30/113__________                          31/113__________                          32/113__________                          33/113__________                          34/113__________                          35/113__________                          36/113__________                          37/113__________                          38/113__________                          39/113__________                          40/113__________                          41/113__________                          42/113__________                          43/113__________                          44/113__________                          45/113__________                          46/113__________                          47/113__________                          48/113__________                          49/113__________                          50/113__________                          51/113__________                          52/113__________                          53/113__________                          54/113__________                          55/113__________                          56/113__________                          57/113__________                          58/113__________                          59/113__________                          60/113__________                          61/113__________                          62/113__________                          63/113__________                          64/113__________                          65/113__________                          66/113__________                          67/113__________                          68/113__________                          69/113__________                          70/113__________                          71/113__________                          72/113__________                          73/113__________                          74/113__________                          75/113__________                          76/113__________                          77/113__________                          78/113__________                          79/113__________                          80/113__________                          81/113__________                          82/113__________                          83/113__________                          84/113__________                          85/113__________                          86/113__________                          87/113__________                          88/113__________                          89/113__________                          90/113__________                          91/113__________                          92/113__________                          93/113__________                          94/113__________                          95/113__________                          96/113__________                          97/113__________                          98/113__________                          99/113__________                          100/113__________                          101/113__________                          102/113__________                          103/113__________                          104/113__________                          105/113__________                          106/113__________                          107/113__________                          108/113__________                          109/113__________                          110/113__________                          111/113__________                          112/113__________                          113/113__________
plot(sbl_2)
#> GH distance not available in this object.


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
#> Predicting...
#> 1/113__________                          2/113__________                          3/113__________                          4/113__________                          5/113__________                          6/113__________                          7/113__________                          8/113__________                          9/113__________                          10/113__________                          11/113__________                          12/113__________                          13/113__________                          14/113__________                          15/113__________                          16/113__________                          17/113__________                          18/113__________                          19/113__________                          20/113__________                          21/113__________                          22/113__________                          23/113__________                          24/113__________                          25/113__________                          26/113__________                          27/113__________                          28/113__________                          29/113__________                          30/113__________                          31/113__________                          32/113__________                          33/113__________                          34/113__________                          35/113__________                          36/113__________                          37/113__________                          38/113__________                          39/113__________                          40/113__________                          41/113__________                          42/113__________                          43/113__________                          44/113__________                          45/113__________                          46/113__________                          47/113__________                          48/113__________                          49/113__________                          50/113__________                          51/113__________                          52/113__________                          53/113__________                          54/113__________                          55/113__________                          56/113__________                          57/113__________                          58/113__________                          59/113__________                          60/113__________                          61/113__________                          62/113__________                          63/113__________                          64/113__________                          65/113__________                          66/113__________                          67/113__________                          68/113__________                          69/113__________                          70/113__________                          71/113__________                          72/113__________                          73/113__________                          74/113__________                          75/113__________                          76/113__________                          77/113__________                          78/113__________                          79/113__________                          80/113__________                          81/113__________                          82/113__________                          83/113__________                          84/113__________                          85/113__________                          86/113__________                          87/113__________                          88/113__________                          89/113__________                          90/113__________                          91/113__________                          92/113__________                          93/113__________                          94/113__________                          95/113__________                          96/113__________                          97/113__________                          98/113__________                          99/113__________                          100/113__________                          101/113__________                          102/113__________                          103/113__________                          104/113__________                          105/113__________                          106/113__________                          107/113__________                          108/113__________                          109/113__________                          110/113__________                          111/113__________                          112/113__________                          113/113__________
plot(local_algo)
#> GH distance not available in this object.


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
#> Predicting...
#> 1/113__________                          2/113__________                          3/113__________                          4/113__________                          5/113__________                          6/113__________                          7/113__________                          8/113__________                          9/113__________                          10/113__________                          11/113__________                          12/113__________                          13/113__________                          14/113__________                          15/113__________                          16/113__________                          17/113__________                          18/113__________                          19/113__________                          20/113__________                          21/113__________                          22/113__________                          23/113__________                          24/113__________                          25/113__________                          26/113__________                          27/113__________                          28/113__________                          29/113__________                          30/113__________                          31/113__________                          32/113__________                          33/113__________                          34/113__________                          35/113__________                          36/113__________                          37/113__________                          38/113__________                          39/113__________                          40/113__________                          41/113__________                          42/113__________                          43/113__________                          44/113__________                          45/113__________                          46/113__________                          47/113__________                          48/113__________                          49/113__________                          50/113__________                          51/113__________                          52/113__________                          53/113__________                          54/113__________                          55/113__________                          56/113__________                          57/113__________                          58/113__________                          59/113__________                          60/113__________                          61/113__________                          62/113__________                          63/113__________                          64/113__________                          65/113__________                          66/113__________                          67/113__________                          68/113__________                          69/113__________                          70/113__________                          71/113__________                          72/113__________                          73/113__________                          74/113__________                          75/113__________                          76/113__________                          77/113__________                          78/113__________                          79/113__________                          80/113__________                          81/113__________                          82/113__________                          83/113__________                          84/113__________                          85/113__________                          86/113__________                          87/113__________                          88/113__________                          89/113__________                          90/113__________                          91/113__________                          92/113__________                          93/113__________                          94/113__________                          95/113__________                          96/113__________                          97/113__________                          98/113__________                          99/113__________                          100/113__________                          101/113__________                          102/113__________                          103/113__________                          104/113__________                          105/113__________                          106/113__________                          107/113__________                          108/113__________                          109/113__________                          110/113__________                          111/113__________                          112/113__________                          113/113__________
plot(local_algo_2)
#> GH distance not available in this object.


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
#> Predicting...
#> 
#> Error in get(x, envir, mode, inherits): argument "method" is missing, with no default

registerDoSEQ()
try(stopCluster(clust))
# }
```
