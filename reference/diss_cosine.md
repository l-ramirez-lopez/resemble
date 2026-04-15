# Cosine dissimilarity method constructor

Creates a configuration object for computing cosine dissimilarity (also
known as spectral angle mapper). Pass the result to
[`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md)
to compute the dissimilarity matrix.

The cosine dissimilarity between two observations \\x_i\\ and \\x_j\\
is:

\$\$c(x_i, x_j) = \cos^{-1} \frac{\sum\_{k=1}^{p} x\_{i,k}\\ x\_{j,k}}
{\sqrt{\sum\_{k=1}^{p} x\_{i,k}^{2}}\\ \sqrt{\sum\_{k=1}^{p}
x\_{j,k}^{2}}}\$\$

where \\p\\ is the number of variables.

## Usage

``` r
diss_cosine(center = TRUE, scale = FALSE)
```

## Arguments

- center:

  Logical. Center the data before computing dissimilarities? Applied
  jointly to `Xr` and `Xu` if both are provided. Default `TRUE`.

- scale:

  Logical. Scale the data before computing dissimilarities? Applied
  jointly to `Xr` and `Xu` if both are provided. Default `FALSE`.

## Value

An object of class `c("diss_cosine", "diss_method")`.

## See also

[`dissimilarity`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md),
[`diss_euclidean`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md),
[`diss_mahalanobis`](https://l-ramirez-lopez.github.io/resemble/reference/diss_mahalanobis.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
m <- diss_cosine()
m <- diss_cosine(center = FALSE)
```
