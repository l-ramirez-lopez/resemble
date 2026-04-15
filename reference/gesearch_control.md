# Control parameters for gesearch

Creates a control object specifying algorithm parameters for
[`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md).

## Usage

``` r
gesearch_control(
  retain_by = c("probability", "proportion"),
  percentile_type = 7L,
  tune = FALSE,
  number = 10L,
  p = 0.75,
  stagnation_limit = 5L,
  allow_parallel = TRUE,
  blas_threads = 1L
)
```

## Arguments

- retain_by:

  A character string specifying how training observations are selected
  at each iteration:

  `"probability"`

  :   (default) Retains observations with errors below a percentile
      estimated from a given probability. More robust to outliers.

  `"proportion"`

  :   Retains a fixed proportion of observations with lowest errors.

- percentile_type:

  An integer between 1 and 9 specifying the quantile algorithm when
  `retain_by = "probability"`. Passed to the `type` argument of
  [`quantile`](https://rdrr.io/r/stats/quantile.html). Default is 7.
  Ignored when `retain_by = "proportion"`.

- tune:

  A logical indicating whether to tune regression parameters (e.g.,
  number of PLS components) via cross-validation at each iteration.
  Increases computation time substantially. Default is `FALSE`.

- number:

  An integer specifying the number of groups for leave-group-out
  cross-validation. Default is 10. This is used for validating the final
  models built the samples found. When `tune = TRUE`, this value is also
  used in the internal CV for tuning regression parameters at each
  iteration.

- p:

  A numeric value in (0, 1) specifying the proportion of observations
  per group in leave-group-out cross-validation. Default is 0.75. When
  `tune = TRUE`, this value is also used in the internal CV for tuning
  regression parameters at each iteration.

- stagnation_limit:

  An integer specifying the maximum number of consecutive iterations
  with no change in gene pool size before early termination. Prevents
  infinite loops when target size cannot be reached. Default is 5.

- allow_parallel:

  A logical indicating whether to enable parallel processing for
  internal resampling and calibration. The parallel backend must be
  registered by the user. Default is `TRUE`.

- blas_threads:

  An integer specifying the number of BLAS threads to use during
  computation. Default is 1, which avoids multi-threaded OpenBLAS
  overhead on Linux. Requires RhpcBLASctl. See Details.

## Value

A list of class `"gesearch_control"` containing the specified
parameters.

## Details

### Retention strategies

When `retain_by = "probability"` (default), observations with errors
below a percentile threshold are retained. The percentile is computed
using [`quantile`](https://rdrr.io/r/stats/quantile.html) with `probs`
set to the `retain` value from
[`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md).
This approach is more robust when outlier observations have extreme
error values.

When `retain_by = "proportion"`, a fixed fraction of observations
(specified by the `retain` argument in
[`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md))
with the lowest associated errors are kept at each iteration.

### Cross-validation for tuning

When `tune = TRUE`, leave-group-out cross-validation is used to select
optimal regression parameters at each iteration. The `number` argument
controls how many CV groups are formed, and `p` controls the proportion
of observations in each group.

### BLAS threading

On Linux systems with multi-threaded OpenBLAS, the default thread count
can cause significant overhead for algorithms that perform many small
matrix operations (like the iterative PLS fits in `gesearch`). Setting
`blas_threads = 1` (the default) eliminates this overhead.

This setting requires the RhpcBLASctl package. If not installed, the
parameter is ignored and a message is displayed. The original thread
count is restored when
[`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md)
completes.

## References

Lobsey, C.R., Viscarra Rossel, R.A., Roudier, P., Hedley, C.B. 2017.
rs-local data-mines information from spectral libraries to improve local
calibrations. European Journal of Soil Science 68:840-852.

## See also

[`gesearch`](https://l-ramirez-lopez.github.io/resemble/reference/gesearch.md),
[`quantile`](https://rdrr.io/r/stats/quantile.html)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# Default parameters (probability-based retention)
gesearch_control()
#> $retain_by
#> [1] "probability"
#> 
#> $percentile_type
#> [1] 7
#> 
#> $tune
#> [1] FALSE
#> 
#> $number
#> [1] 10
#> 
#> $p
#> [1] 0.75
#> 
#> $stagnation_limit
#> [1] 5
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "gesearch_control"

# Proportion-based retention
gesearch_control(retain_by = "proportion")
#> $retain_by
#> [1] "proportion"
#> 
#> $percentile_type
#> NULL
#> 
#> $tune
#> [1] FALSE
#> 
#> $number
#> [1] 10
#> 
#> $p
#> [1] 0.75
#> 
#> $stagnation_limit
#> [1] 5
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "gesearch_control"

# Enable parameter tuning with custom CV settings
gesearch_control(tune = TRUE, number = 5, p = 0.8)
#> $retain_by
#> [1] "probability"
#> 
#> $percentile_type
#> [1] 7
#> 
#> $tune
#> [1] TRUE
#> 
#> $number
#> [1] 5
#> 
#> $p
#> [1] 0.8
#> 
#> $stagnation_limit
#> [1] 5
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "gesearch_control"
```
