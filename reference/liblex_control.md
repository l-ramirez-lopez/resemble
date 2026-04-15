# Control parameters for liblex

Specifies control parameters for the
[`liblex`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md)
function, including output options, validation settings, tuning
behavior, and parallel execution.

## Usage

``` r
liblex_control(
  return_dissimilarity = FALSE,
  mode = c("build", "validate"),
  tune = FALSE,
  metric = c("rmse", "r2"),
  chunk_size = 1L,
  allow_parallel = TRUE,
  blas_threads = 1L
)
```

## Arguments

- return_dissimilarity:

  A logical indicating whether the dissimilarity matrix should be
  returned in the output. Default is `FALSE`. Setting to `TRUE` can be
  useful for diagnostics but increases memory usage for large libraries.

- mode:

  A character string specifying the operation mode:

  `"build"`

  :   (default) Builds the library of local models. If `tune = TRUE`,
      validation is performed first to find optimal parameters, then the
      library is built using those parameters. If `tune = FALSE`, the
      library is built directly using the parameters provided to
      [`liblex`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md).

  `"validate"`

  :   Performs validation only without building the library. Useful for
      parameter exploration or when testing different configurations
      before committing to the full library build.

- tune:

  A logical indicating whether to optimize parameters via
  nearest-neighbor validation. Default is `FALSE`. When `TRUE`, the
  function evaluates all combinations of `k` values and PLS component
  ranges, selecting the combination that minimizes RMSE (or maximizes
  R², depending on `metric`). See Details.

- metric:

  A character string specifying the performance metric used for
  parameter selection when `tune = TRUE`. Options are:

  `"rmse"`

  :   (default) Root mean squared error (minimized).

  `"r2"`

  :   Coefficient of determination (maximized).

  Ignored when `tune = FALSE`.

- chunk_size:

  An integer specifying the number of local models to process per
  parallel task. Default is `1L`. Increasing this value reduces parallel
  overhead but may cause load imbalance. For large libraries, values
  between 10 and 50 often provide a good trade-off between overhead and
  efficiency.

- allow_parallel:

  A logical indicating whether parallel execution is permitted. Default
  is `TRUE`. Parallelization is applied to the model fitting loop using
  the foreach package. Requires a registered parallel backend (e.g., via
  doParallel).

- blas_threads:

  An integer specifying the number of threads for BLAS operations.
  Default is `1L`, which avoids thread contention when using parallel
  processing. Requires the RhpcBLASctl package to take effect. On Linux
  systems with multi-threaded BLAS (e.g., OpenBLAS), setting this to 1
  can substantially improve performance when `allow_parallel = TRUE`.

## Value

A list of class `liblex_control` containing the validated control
parameters.

## Details

### Nearest-neighbor validation

When `tune = TRUE` or `mode = "validate"`, the function performs
nearest-neighbor validation (NNv) to assess model performance. For each
observation in the reference set (or anchor set, if specified), the
procedure:

1.  Identifies the k nearest neighbors of the target observation.

2.  Excludes the target observation (and any observations in the same
    group, if `group` is specified) from the neighbor set.

3.  Fits a local model using the remaining neighbors.

4.  Predicts the response value for the excluded target observation.

5.  Computes prediction errors across all observations.

This leave-one-out style validation provides an estimate of prediction
performance without requiring a separate test set. When `tune = TRUE`,
the parameter combination (number of neighbors and PLS component range)
yielding the best performance according to `metric` is selected for
building the final library.

### Mode and tune combinations

|              |         |                                                                |
|--------------|---------|----------------------------------------------------------------|
| `mode`       | `tune`  | Behavior                                                       |
| `"build"`    | `FALSE` | Build library using parameters as provided                     |
| `"build"`    | `TRUE`  | Validate, find optimal parameters, build library               |
| `"validate"` | `FALSE` | Validate only, report performance statistics                   |
| `"validate"` | `TRUE`  | Validate, report statistics with optimal parameters identified |

### Parallel chunk size

The `chunk_size` parameter controls granularity of parallel work
distribution. When `allow_parallel = TRUE` and a parallel backend is
registered:

- `chunk_size = 1`: Each local model is a separate parallel task.
  Maximum parallelism but higher scheduling overhead.

- `chunk_size > 1`: Multiple models are processed sequentially within
  each parallel task. Reduces overhead and improves memory locality, but
  may cause load imbalance if the number of models is not evenly
  divisible.

When `allow_parallel = FALSE`, `chunk_size` has no effect.

## See also

[`liblex`](https://l-ramirez-lopez.github.io/resemble/reference/liblex.md),
[`mbl_control`](https://l-ramirez-lopez.github.io/resemble/reference/mbl_control.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# Default settings: build library without tuning
liblex_control()
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $mode
#> [1] "build"
#> 
#> $tune
#> [1] FALSE
#> 
#> $metric
#> [1] "rmse"
#> 
#> $chunk_size
#> [1] 1
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "liblex_control"

# Tune parameters before building
liblex_control(tune = TRUE)
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $mode
#> [1] "build"
#> 
#> $tune
#> [1] TRUE
#> 
#> $metric
#> [1] "rmse"
#> 
#> $chunk_size
#> [1] 1
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "liblex_control"

# Validate only for parameter exploration
liblex_control(mode = "validate", tune = TRUE)
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $mode
#> [1] "validate"
#> 
#> $tune
#> [1] TRUE
#> 
#> $metric
#> [1] "rmse"
#> 
#> $chunk_size
#> [1] 1
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "liblex_control"

# Include dissimilarity matrix in output
liblex_control(return_dissimilarity = TRUE)
#> $return_dissimilarity
#> [1] TRUE
#> 
#> $mode
#> [1] "build"
#> 
#> $tune
#> [1] FALSE
#> 
#> $metric
#> [1] "rmse"
#> 
#> $chunk_size
#> [1] 1
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "liblex_control"

# Larger chunks for reduced parallel overhead
liblex_control(chunk_size = 20L)
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $mode
#> [1] "build"
#> 
#> $tune
#> [1] FALSE
#> 
#> $metric
#> [1] "rmse"
#> 
#> $chunk_size
#> [1] 20
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "liblex_control"

# Parallel settings for Linux with OpenBLAS
liblex_control(allow_parallel = TRUE, blas_threads = 1L)
#> $return_dissimilarity
#> [1] FALSE
#> 
#> $mode
#> [1] "build"
#> 
#> $tune
#> [1] FALSE
#> 
#> $metric
#> [1] "rmse"
#> 
#> $chunk_size
#> [1] 1
#> 
#> $allow_parallel
#> [1] TRUE
#> 
#> $blas_threads
#> [1] 1
#> 
#> attr(,"class")
#> [1] "liblex_control"
```
