# Mahalanobis dissimilarity method constructor

Creates a configuration object for computing Mahalanobis dissimilarity.
Pass the result to
[`dissimilarity()`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md)
to compute the dissimilarity matrix.

The Mahalanobis distance is computed by first transforming the data into
Mahalanobis space via a factorization of the inverse covariance matrix
\\M^{-1} = W^{T}W\\ (using SVD), then applying Euclidean distance in
that transformed space:

\\d(x_i, x_j) = \sqrt{\frac{1}{p}(x_i - x_j)M^{-1}(x_i - x_j)^T}\\

## Usage

``` r
diss_mahalanobis(center = TRUE, scale = FALSE)
```

## Arguments

- center:

  Logical. Center the data before computing distances? Applied jointly
  to `Xr` and `Xu` if both are provided. Default `TRUE`.

- scale:

  Logical. Scale the data before computing distances? Applied jointly to
  `Xr` and `Xu` if both are provided. Default `FALSE`.

## Value

An object of class `c("diss_mahalanobis", "diss_method")`.

## Important limitations

The covariance matrix will be singular — and the distance therefore
uncomputable — when the number of observations is smaller than the
number of variables, or when variables are perfectly collinear. This is
common with raw spectral data; consider using
[`diss_euclidean()`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md)
on PCA scores instead.

## See also

[`dissimilarity`](https://l-ramirez-lopez.github.io/resemble/reference/dissimilarity.md),
[`diss_euclidean`](https://l-ramirez-lopez.github.io/resemble/reference/diss_euclidean.md),
[`diss_cosine`](https://l-ramirez-lopez.github.io/resemble/reference/diss_cosine.md)

## Author

[Leonardo Ramirez-Lopez](https://orcid.org/0000-0002-5369-5120)

## Examples

``` r
m <- diss_mahalanobis()
m <- diss_mahalanobis(center = TRUE, scale = TRUE)
```
