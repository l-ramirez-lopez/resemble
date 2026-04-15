# Control parameters for memory-based learning

This function controls various aspects of the memory-based learning
process in the
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
function.

## Usage

``` r
mbl_control(
  return_dissimilarity = FALSE,
  validation_type = "NNv",
  tune_locally = TRUE,
  number = 10,
  p = 0.75,
  range_prediction_limits = TRUE,
  allow_parallel = TRUE,
  blas_threads = 1L
)
```

## Arguments

- return_dissimilarity:

  Logical indicating whether to return the dissimilarity matrix between
  `Xr` and `Xu`. Default is `FALSE`.

- validation_type:

  Character vector specifying validation method(s):

  - `"NNv"`: Leave-nearest-neighbor-out cross-validation (default,
    faster)

  - `"local_cv"`: Local leave-group-out cross-validation

  - `"none"`: No validation

  Multiple methods can be specified (e.g., `c("NNv", "local_cv")`).
  Default is `"NNv"`.

- tune_locally:

  Logical indicating whether to tune PLS components locally when
  `validation_type = "local_cv"` and using
  [`fit_pls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md)
  or
  [`fit_wapls`](https://l-ramirez-lopez.github.io/resemble/reference/fit_methods.md).
  Default is `TRUE`.

- number:

  Integer specifying the number of sampling iterations for `"local_cv"`
  validation. Default is `10`.

- p:

  Numeric value between 0 and 1 indicating the proportion of
  observations retained at each `"local_cv"` iteration. Default is
  `0.75`.

- range_prediction_limits:

  Logical indicating whether predictions should be constrained to the
  range of response values in each neighborhood. Default is `TRUE`.

- allow_parallel:

  Logical indicating whether parallel execution is allowed via the
  foreach package. Default is `TRUE`.

- blas_threads:

  Integer specifying the number of BLAS threads to use during
  [`mbl()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
  execution. Default is `1L`, which avoids thread overhead from repeated
  small matrix operations. Requires the RhpcBLASctl package to take
  effect. The original thread count is restored after
  [`mbl()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
  completes. See Details.

## Value

A list of class `mbl_control` with the specified control parameters.

## Details

### Validation methods

**Leave-nearest-neighbor-out cross-validation (`"NNv"`):** For each
target observation, the nearest neighbor is excluded from the local
model, which then predicts that neighbor's value. This is faster than
`"local_cv"`. If the nearest neighbor belongs to a group (specified via
the `group` argument in
[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)),
all group members are excluded.

**Local leave-group-out cross-validation (`"local_cv"`):** The
neighborhood is partitioned into subsets via stratified random sampling.
Each subset serves as validation data while the remainder fits the
model. This repeats `number` times, with `p` controlling the training
proportion. The final error is the average local RMSE.

### BLAS threading

On Linux systems with multi-threaded OpenBLAS, the default thread count
can cause significant overhead for algorithms like
[`mbl()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
that perform many small matrix operations. Setting `blas_threads = 1`
(the default) eliminates this overhead.

This setting requires the RhpcBLASctl package. If not installed, the
parameter is ignored and a message is displayed. The original thread
count is restored when
[`mbl()`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md)
completes.

Windows systems typically use single-threaded BLAS by default, so this
setting has no effect there.

## References

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte,
J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local
approach for modeling soil vis-NIR spectra of complex data sets.
Geoderma 195-196:268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R.,
Dematte, J.A.M., Scholten, T. 2013b. Distance and similarity-search
metrics for use with soil vis-NIR spectra. Geoderma 199:43-53.

## See also

[`mbl`](https://l-ramirez-lopez.github.io/resemble/reference/mbl.md),
[`neighbors_k`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md),
[`neighbors_diss`](https://l-ramirez-lopez.github.io/resemble/reference/neighbors.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120) and
Antoine Stevens

## Examples

``` r
# Default control parameters (NNv validation)
mbl_control()
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $validation_type
#> [1] "NNv"
#> 
#> $tune_locally
#> [1] TRUE
#> 
#> $number
#> [1] 10
#> 
#> $p
#> [1] 0.75
#> 
#> $range_prediction_limits
#> [1] TRUE
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "mbl_control" "list"       

# Both validation methods
mbl_control(validation_type = c("NNv", "local_cv"))
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $validation_type
#> [1] "NNv"      "local_cv"
#> 
#> $tune_locally
#> [1] TRUE
#> 
#> $number
#> [1] 10
#> 
#> $p
#> [1] 0.75
#> 
#> $range_prediction_limits
#> [1] TRUE
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "mbl_control" "list"       

# No validation
mbl_control(validation_type = "none")
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $validation_type
#> [1] "none"
#> 
#> $tune_locally
#> [1] TRUE
#> 
#> $number
#> [1] 10
#> 
#> $p
#> [1] 0.75
#> 
#> $range_prediction_limits
#> [1] TRUE
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "mbl_control" "list"       

# NNv validation only, no parallel
mbl_control(validation_type = "NNv", allow_parallel = FALSE)
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $validation_type
#> [1] "NNv"
#> 
#> $tune_locally
#> [1] TRUE
#> 
#> $number
#> [1] 10
#> 
#> $p
#> [1] 0.75
#> 
#> $range_prediction_limits
#> [1] TRUE
#> 
#> $allow_parallel
#> [1] FALSE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "mbl_control" "list"       

# Allow more BLAS threads (if needed for other computations)
mbl_control(blas_threads = 4)
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $validation_type
#> [1] "NNv"
#> 
#> $tune_locally
#> [1] TRUE
#> 
#> $number
#> [1] 10
#> 
#> $p
#> [1] 0.75
#> 
#> $range_prediction_limits
#> [1] TRUE
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 4
#> 
#> attr(,"class")
#> [1] "mbl_control" "list"       
```
