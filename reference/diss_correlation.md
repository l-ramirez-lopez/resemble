# Correlation dissimilarity method constructor

Creates a configuration object that fully specifies a correlation (or
moving correlation) dissimilarity method. Pass the result to
[`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md)
to compute the dissimilarity matrix.

## Usage

``` r
diss_correlation(ws = NULL, center = TRUE, scale = FALSE)
```

## Arguments

- ws:

  Either `NULL` (default) or an odd integer greater than 2. When `NULL`,
  standard Pearson correlation dissimilarity is used. When an odd
  integer is provided, a moving (rolling) correlation dissimilarity is
  computed using a window of that size.

- center:

  Logical. Should the data be mean-centered before computing
  dissimilarities? Centering is applied jointly to `Xr` and `Xu` (if
  provided) based on their combined column means. Default is `TRUE`.

- scale:

  Logical. Should the data be scaled (divided by column standard
  deviations) before computing dissimilarities? Scaling is applied
  jointly to `Xr` and `Xu` (if provided). Default is `FALSE`.

## Value

An object of class `c("diss_correlation", "diss_method")` — a list
holding the validated method parameters. Intended to be passed to
[`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md),
not used directly.

## Details

The correlation dissimilarity between two observations \\x_i\\ and
\\x_j\\ is:

\$\$d(x_i, x_j) = \frac{1}{2}(1 - \rho(x_i, x_j))\$\$

where \\\rho\\ is the Pearson correlation coefficient. This is used when
`ws = NULL`.

When `ws` is specified, the moving correlation dissimilarity is:

\$\$d(x_i, x_j; ws) = \frac{1}{2\\ws} \sum\_{k=1}^{p - ws} \bigl(1 -
\rho(x\_{i,(k:k+ws)},\\ x\_{j,(k:k+ws)})\bigr)\$\$

where \\ws\\ is the window size and \\p\\ is the number of variables.

## Parallel execution

The underlying C++ implementation uses OpenMP for parallel computation.
Thread count is controlled by the `OMP_NUM_THREADS` environment
variable. To limit threads (e.g., when calling from within a parallel
backend):

    Sys.setenv(OMP_NUM_THREADS = 1)

## See also

[`dissimilarity`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md),
[`diss_euclidean`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md),
[`diss_mahalanobis`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md),
[`diss_cosine`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
# Standard correlation dissimilarity
m <- diss_correlation()

# Moving correlation with window size 41
m <- diss_correlation(ws = 41)

# Without centering
m <- diss_correlation(center = FALSE)
```
