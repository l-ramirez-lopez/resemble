# Euclidean dissimilarity method constructor

Creates a configuration object for computing Euclidean dissimilarity.
Pass the result to
[`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md)
to compute the dissimilarity matrix.

The scaled Euclidean dissimilarity between two observations \\x_i\\ and
\\x_j\\ is:

\$\$d(x_i, x_j) = \sqrt{\frac{1}{p} \sum\_{k=1}^{p}(x\_{i,k} -
x\_{j,k})^2}\$\$

where \\p\\ is the number of variables. Results are equivalent to
[`stats::dist()`](https://rdrr.io/r/stats/dist.html) but scaled by
\\1/p\\.

## Usage

``` r
diss_euclidean(center = TRUE, scale = FALSE)
```

## Arguments

- center:

  Logical. Center the data before computing distances? Applied jointly
  to `Xr` and `Xu` if both are provided. Default `TRUE`.

- scale:

  Logical. Scale the data before computing distances? Applied jointly to
  `Xr` and `Xu` if both are provided. Default `FALSE`.

## Value

An object of class `c("diss_euclidean", "diss_method")`.

## See also

[`dissimilarity`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md),
[`diss_mahalanobis`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md),
[`diss_cosine`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
m <- diss_euclidean()
m <- diss_euclidean(center = FALSE, scale = TRUE)
```
