context("test-euclidean")

test_that("dist() equivalence to resemble euclidean", {

  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  tol <- 1e-3
  dist_m <- as.matrix(dist(Xr, method = "euclidean"))
  diss_m <- dissimilarity(
    Xr = Xr, diss_method = diss_euclidean(center = FALSE, scale = FALSE)
  )$dissimilarity
  diss_m <- sqrt(diss_m^2 * ncol(Xr))
  expect_lt(sum(abs(diss_m - dist_m)), tol)
})

