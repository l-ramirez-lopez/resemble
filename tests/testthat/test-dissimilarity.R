context("test-dissimilarity")

# tests/testthat/test-dissimilarity.R
library(foreach)
library(RhpcBLASctl)
registerDoSEQ()
# =============================================================================
# Setup helper
# =============================================================================

.setup_nirsoil_data <- function(n_xr = 40, n_xu = 20) {
  data("NIRsoil", package = "prospectr")
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ][seq_len(n_xu), ]
  Xr <- Xr[!is.na(Yr), ][seq_len(n_xr), ]
  Yu <- Yu[!is.na(Yu)][seq_len(n_xu)]
  Yr <- Yr[!is.na(Yr)][seq_len(n_xr)]
  
  list(Xr = Xr, Xu = Xu, Yr = Yr, Yu = Yu)
}



# =============================================================================
# diss_* constructor tests
# =============================================================================

test_that("diss_correlation constructor works", {
  x <- diss_correlation(ws = 11)
  expect_s3_class(x, "diss_method")
})

# =============================================================================
# Tests for diss_euclidean, diss_mahalanobis, diss_cosine constructors
# =============================================================================


# -----------------------------------------------------------------------------
# diss_euclidean constructor tests
# -----------------------------------------------------------------------------

test_that("diss_euclidean returns correct class", {
  m <- diss_euclidean()
  expect_s3_class(m, "diss_euclidean")
  expect_s3_class(m, "diss_method")
})

test_that("diss_euclidean stores parameters correctly", {
  m1 <- diss_euclidean()
  expect_true(m1$center)
  expect_false(m1$scale)
  expect_equal(m1$method, "euclidean")
  
  m2 <- diss_euclidean(center = FALSE, scale = TRUE)
  expect_false(m2$center)
  expect_true(m2$scale)
})

test_that("diss_euclidean computes correct dissimilarities", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  
  # Test without centering/scaling - compare to stats::dist
  result <- dissimilarity(Xr, diss_method = diss_euclidean(center = FALSE, scale = FALSE))
  
  expect_true(is.matrix(result$dissimilarity))
  expect_equal(nrow(result$dissimilarity), nrow(Xr))
  expect_equal(ncol(result$dissimilarity), nrow(Xr))
  
  # Diagonal should be zero (self-dissimilarity)
  expect_true(all(abs(diag(result$dissimilarity)) < 1e-6))
  
  # Should be symmetric
  expect_equal(result$dissimilarity, t(result$dissimilarity), tolerance = 1e-10)
})

test_that("diss_euclidean with Xu computes cross-dissimilarity", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  Xu <- NIRsoil$spc[21:30, ]
  
  
  result <- dissimilarity(Xr, Xu, diss_method = diss_euclidean())
  
  expect_equal(nrow(result$dissimilarity), nrow(Xr))
  expect_equal(ncol(result$dissimilarity), nrow(Xu))
})

test_that("diss_euclidean centering affects results", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  
  result_centered <- dissimilarity(Xr, diss_method = diss_euclidean(center = TRUE))
  
  result_uncentered <- dissimilarity(Xr, diss_method = diss_euclidean(center = FALSE))
  
  # Results should differ when centering is applied
  expect_false(all(result_centered$dissimilarity == result_uncentered$dissimilarity))
})

test_that("diss_euclidean scaling affects results", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  
  result_scaled <- dissimilarity(Xr, diss_method = diss_euclidean(scale = TRUE))
  result_unscaled <- dissimilarity(Xr, diss_method = diss_euclidean(scale = FALSE))
  
  # Results should differ when scaling is applied
  expect_false(all(result_scaled$dissimilarity == result_unscaled$dissimilarity))
})

# -----------------------------------------------------------------------------
# diss_mahalanobis constructor tests
# -----------------------------------------------------------------------------

test_that("diss_mahalanobis returns correct class", {
  m <- diss_mahalanobis()
  expect_s3_class(m, "diss_mahalanobis")
  expect_s3_class(m, "diss_method")
})

test_that("diss_mahalanobis stores parameters correctly", {
  m1 <- diss_mahalanobis()
  expect_true(m1$center)
  expect_false(m1$scale)
  expect_equal(m1$method, "mahalanobis")
  
  m2 <- diss_mahalanobis(center = FALSE, scale = TRUE)
  expect_false(m2$center)
  expect_true(m2$scale)
})

test_that("diss_mahalanobis computes dissimilarities when n > p", {
  # Create data where n > p (more observations than variables)
  skip_on_cran()
  set.seed(123)
  Xr <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  
  result <- dissimilarity(Xr, diss_method = diss_mahalanobis())
  
  expect_true(is.matrix(result$dissimilarity))
  expect_equal(nrow(result$dissimilarity), nrow(Xr))
  expect_equal(ncol(result$dissimilarity), nrow(Xr))
  
  # Diagonal should be zero
  expect_true(all(abs(diag(result$dissimilarity)) < 1e-6))
  
  # Should be symmetric
  expect_equal(result$dissimilarity, t(result$dissimilarity), tolerance = 1e-10)
})

test_that("diss_mahalanobis errors when covariance is singular", {
  # Create data where n < p (fewer observations than variables) - singular covariance
  set.seed(123)
  Xr <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
  
  expect_error(
    dissimilarity(Xr, diss_method = diss_mahalanobis()),
    "For computing the Mahalanobis distance, the total number of observations"
  )
})

test_that("diss_mahalanobis with Xu computes cross-dissimilarity", {
  skip_on_cran()
  set.seed(123)
  Xr <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  Xu <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
  
  result <- dissimilarity(Xr, Xu, diss_method = diss_mahalanobis())
  
  expect_equal(nrow(result$dissimilarity), nrow(Xr))
  expect_equal(ncol(result$dissimilarity), nrow(Xu))
})

test_that("diss_mahalanobis differs from diss_euclidean", {
  set.seed(123)
  Xr <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  
  result_mahal <- dissimilarity(Xr, diss_method = diss_mahalanobis())
  result_euclid <- dissimilarity(Xr, diss_method = diss_euclidean())
  
  # Mahalanobis and Euclidean should give different results
  expect_false(all(result_mahal$dissimilarity == result_euclid$dissimilarity))
})

# -----------------------------------------------------------------------------
# diss_cosine constructor tests
# -----------------------------------------------------------------------------

test_that("diss_cosine returns correct class", {
  m <- diss_cosine()
  expect_s3_class(m, "diss_cosine")
  expect_s3_class(m, "diss_method")
})

test_that("diss_cosine stores parameters correctly", {
  m1 <- diss_cosine()
  expect_true(m1$center)
  expect_false(m1$scale)
  expect_equal(m1$method, "cosine")
  
  m2 <- diss_cosine(center = FALSE, scale = TRUE)
  expect_false(m2$center)
  expect_true(m2$scale)
})

test_that("diss_cosine computes correct dissimilarities", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  
  result <- dissimilarity(Xr, diss_method = diss_cosine())
  
  expect_true(is.matrix(result$dissimilarity))
  expect_equal(nrow(result$dissimilarity), nrow(Xr))
  expect_equal(ncol(result$dissimilarity), nrow(Xr))
  
  # Diagonal should be zero (self-dissimilarity)
  expect_true(all(abs(diag(result$dissimilarity)) < 1e-6))
  
  # Should be symmetric
  expect_equal(result$dissimilarity, t(result$dissimilarity), tolerance = 1e-10)
})

test_that("diss_cosine values are bounded", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  
  result <- dissimilarity(Xr, diss_method = diss_cosine(center = FALSE))
  
  # Cosine dissimilarity (angle) should be between 0 and pi
  expect_true(all(result$dissimilarity >= 0))
  expect_true(all(result$dissimilarity <= pi + 1e-10))
})

test_that("diss_cosine with Xu computes cross-dissimilarity", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  Xu <- NIRsoil$spc[21:30, ]
  
  result <- dissimilarity(Xr, Xu, diss_method = diss_cosine())
  
  expect_equal(nrow(result$dissimilarity), nrow(Xr))
  expect_equal(ncol(result$dissimilarity), nrow(Xu))
})

test_that("diss_cosine differs from diss_euclidean", {
  skip_if_not_installed("prospectr")
  data("NIRsoil", package = "prospectr")
  Xr <- NIRsoil$spc[1:20, ]
  
  result_cosine <- dissimilarity(Xr, diss_method = diss_cosine())
  result_euclid <- dissimilarity(Xr, diss_method = diss_euclidean())
  
  # Cosine and Euclidean should give different results
  expect_false(all(result_cosine$dissimilarity == result_euclid$dissimilarity))
})

test_that("diss_cosine handles identical observations", {
  # Two identical rows should have zero dissimilarity
  Xr <- matrix(c(1, 2, 3, 1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
  
  result <- dissimilarity(Xr, diss_method = diss_cosine(center = FALSE))
  
  # First two rows are identical
  
  expect_equal(result$dissimilarity[1, 2], 0, tolerance = 1e-10)
  expect_equal(result$dissimilarity[2, 1], 0, tolerance = 1e-10)
})

# -----------------------------------------------------------------------------
# Print method tests
# -----------------------------------------------------------------------------

test_that("print.diss_method works for all constructors", {
  expect_output(print(diss_euclidean()), "euclidean")
  expect_output(print(diss_mahalanobis()), "mahalanobis")
  expect_output(print(diss_cosine()), "cosine")
})

# -----------------------------------------------------------------------------
# Edge cases
# -----------------------------------------------------------------------------

test_that("not possible to have a single observation", {
  Xr <- matrix(1:10, nrow = 1)
  
  expect_error(dissimilarity(Xr, diss_method = diss_euclidean()))
})

test_that("constructors handle two observations", {
  Xr <- matrix(1:20, nrow = 2)
  
  result_euclid <- dissimilarity(Xr, diss_method = diss_euclidean())
  expect_equal(dim(result_euclid$dissimilarity), c(2, 2))
  expect_equal(result_euclid$dissimilarity[1, 2], result_euclid$dissimilarity[2, 1])
  
  result_cosine <- dissimilarity(Xr, diss_method = diss_cosine())
  expect_equal(dim(result_cosine$dissimilarity), c(2, 2))
  expect_equal(result_cosine$dissimilarity[1, 2], result_cosine$dissimilarity[2, 1])
})



test_that("diss_pca constructor works", {
  m <- diss_pca()
  expect_s3_class(m, c("diss_pca", "diss_method"))
  expect_s3_class(m$ncomp, "ncomp_by_var")  # default
  expect_equal(m$method, "pca")
  expect_true(m$center)
  expect_false(m$scale)
  expect_false(m$return_projection)
  
  m2 <- diss_pca(
    ncomp = 10, 
    method = "pca_nipals", 
    center = FALSE, 
    scale = TRUE,
    return_projection = TRUE
  )
  expect_s3_class(m2$ncomp, "ncomp_fixed")
  expect_equal(m2$method, "pca_nipals")
  expect_false(m2$center)
  expect_true(m2$scale)
  expect_true(m2$return_projection)
})


test_that("diss_pls constructor works", {
  m <- diss_pls()
  expect_s3_class(m, c("diss_pls", "diss_method"))
  expect_s3_class(m$ncomp, "ncomp_by_opc")  # default
  expect_equal(m$method, "pls")
  expect_true(m$center)  # PLS always centers
  expect_false(m$scale)
  
  m2 <- diss_pls(ncomp = 15, method = "mpls", scale = TRUE)
  expect_s3_class(m2$ncomp, "ncomp_fixed")
  expect_equal(m2$method, "mpls")
  expect_true(m2$scale)
})


test_that("diss_euclidean constructor works", {
  m <- diss_euclidean()
  expect_s3_class(m, c("diss_euclidean", "diss_method"))
  expect_true(m$center)
  expect_false(m$scale)
  
  m2 <- diss_euclidean(center = FALSE, scale = TRUE)
  expect_false(m2$center)
  expect_true(m2$scale)
})


test_that("diss_mahalanobis constructor works", {
  m <- diss_mahalanobis()
  expect_s3_class(m, c("diss_mahalanobis", "diss_method"))
  expect_true(m$center)
  expect_false(m$scale)
})


test_that("diss_cosine constructor works", {
  m <- diss_cosine()
  expect_s3_class(m, c("diss_cosine", "diss_method"))
  expect_true(m$center)
  expect_false(m$scale)
})


test_that("diss_correlation constructor works", {
  m <- diss_correlation()
  expect_s3_class(m, c("diss_correlation", "diss_method"))
  expect_null(m$ws)
  expect_true(m$center)
  expect_false(m$scale)
  
  m2 <- diss_correlation(ws = 41)
  expect_equal(m2$ws, 41L)
})


test_that("diss_correlation validates ws", {
  expect_error(diss_correlation(ws = 2), "greater than 2")
  expect_error(diss_correlation(ws = 4), "odd integer")
  expect_error(diss_correlation(ws = "a"), "single integer")
})


# =============================================================================
# diss_pca dissimilarity tests
# =============================================================================

test_that("diss_pca with ncomp_by_var works", {
  skip_if_not_installed("prospectr")
  skip_on_cran()
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_pca(
      ncomp = ncomp_by_var(0.01),
      center = TRUE, 
      scale = TRUE,
      return_projection = TRUE
    )
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
  expect_true("projection" %in% names(result))
  expect_true("ncomp" %in% names(result))
  expect_equal(nrow(result$dissimilarity), nrow(d$Xr))
  expect_equal(ncol(result$dissimilarity), nrow(d$Xu))
})


test_that("diss_pca with ncomp_by_opc works", {
  skip_if_not_installed("prospectr")
  skip_on_cran()
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    diss_method = diss_pca(
      ncomp = ncomp_by_opc(max_ncomp = 15),
      center = TRUE, 
      scale = TRUE,
      return_projection = TRUE
    )
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
  expect_true("projection" %in% names(result))
  expect_true(result$ncomp <= 15)
})


test_that("diss_pca with ncomp_by_opc requires Yr", {
  skip_if_not_installed("prospectr")
  skip_on_cran()
  
  d <- .setup_nirsoil_data()
  
  expect_error(
    dissimilarity(
      Xr = d$Xr, 
      Xu = d$Xu,
      diss_method = diss_pca(ncomp = ncomp_by_opc())
    ),
    "'Yr' is required"
  )
})


test_that("diss_pca with ncomp_by_cumvar works", {
  skip_if_not_installed("prospectr")
  skip_on_cran()
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_pca(
      ncomp = ncomp_by_cumvar(0.99),
      return_projection = TRUE
    )
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true(result$ncomp >= 1)
})


test_that("diss_pca with fixed ncomp works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  # Integer form
  result1 <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10)
  )
  
  # Explicit ncomp_fixed form
  result2 <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_pca(ncomp = ncomp_fixed(10))
  )
  
  expect_equal(result1$ncomp, 10)
  expect_equal(result2$ncomp, 10)
  expect_equal(result1$dissimilarity, result2$dissimilarity)
})


test_that("diss_pca with pca_nipals method works", {
  skip_if_not_installed("prospectr")
  skip_on_cran()
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_pca(
      ncomp = 10,
      method = "pca_nipals",
      return_projection = TRUE
    )
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_equal(result$ncomp, 10)
})


test_that("diss_pca self-dissimilarity is symmetric with zero diagonal", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 20, n_xu = 0)
  
  result <- dissimilarity(
    Xr = d$Xr,
    diss_method = diss_pca(ncomp = 5)
  )
  
  expect_equal(result$dissimilarity, t(result$dissimilarity))
  expect_true(all(abs(diag(result$dissimilarity)) < 1e-6))
})


# =============================================================================
# diss_pls dissimilarity tests
# =============================================================================

test_that("diss_pls with ncomp_by_opc works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  Xr_snv <- prospectr::standardNormalVariate(d$Xr)
  Xu_snv <- prospectr::standardNormalVariate(d$Xu)
  
  result <- dissimilarity(
    Xr = Xr_snv, 
    Xu = Xu_snv,
    Yr = d$Yr,
    diss_method = diss_pls(
      ncomp = ncomp_by_opc(max_ncomp = 15),
      return_projection = TRUE
    )
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
  expect_true("projection" %in% names(result))
  expect_true(result$ncomp <= 15)
})


test_that("diss_pls with ncomp_by_var works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    diss_method = diss_pls(
      ncomp = ncomp_by_var(0.01),
      return_projection = TRUE
    )
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true(result$ncomp >= 1)
})


test_that("diss_pls with fixed ncomp works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    diss_method = diss_pls(ncomp = 15)
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_equal(result$ncomp, 15)
})


test_that("diss_pls with mpls method works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    diss_method = diss_pls(
      ncomp = 10,
      method = "mpls"
    )
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_equal(result$ncomp, 10)
})


test_that("diss_pls requires Yr", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  expect_error(
    dissimilarity(
      Xr = d$Xr, 
      Xu = d$Xu,
      diss_method = diss_pls()
    ),
    "'Yr' is required"
  )
})


# =============================================================================
# diss_euclidean dissimilarity tests
# =============================================================================

test_that("diss_euclidean works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_euclidean(center = TRUE, scale = TRUE)
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
  expect_equal(nrow(result$dissimilarity), nrow(d$Xr))
  expect_equal(ncol(result$dissimilarity), nrow(d$Xu))
})


test_that("diss_euclidean matches stats::dist", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 10, n_xu = 0)
  
  # Compute with resemble (no centering/scaling for direct comparison)
  result <- dissimilarity(
    Xr = d$Xr,
    diss_method = diss_euclidean(center = FALSE, scale = FALSE)
  )
  
  # resemble scales by 1/sqrt(p), so convert back
  resemble_diss <- result$dissimilarity
  resemble_diss_unscaled <- ((resemble_diss^2) * ncol(d$Xr))^0.5
  
  # Compare with base R
  base_diss <- as.matrix(dist(d$Xr))
  
  expect_lt(sum(abs(resemble_diss_unscaled - base_diss)), 1e-5)
})


test_that("diss_euclidean self-dissimilarity is symmetric", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 20, n_xu = 0)
  
  result <- dissimilarity(
    Xr = d$Xr,
    diss_method = diss_euclidean()
  )
  
  expect_equal(result$dissimilarity, t(result$dissimilarity))
  expect_true(all(abs(diag(result$dissimilarity)) < 1e-6))
})


# =============================================================================
# diss_correlation dissimilarity tests
# =============================================================================

test_that("diss_correlation (full) works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_correlation(ws = NULL, center = TRUE, scale = FALSE)
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
  expect_equal(nrow(result$dissimilarity), nrow(d$Xr))
  expect_equal(ncol(result$dissimilarity), nrow(d$Xu))
})


test_that("diss_correlation (moving window) works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_correlation(ws = 11)
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
})


test_that("diss_correlation vs one computed with cor()", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_correlation(center = FALSE)
  )
  
  result_base_cor <- 0.5 * (1 - cor(t(d$Xr), t(d$Xu)))
  
  # Results should be similar but not identical 
  expect_lt(
    max(abs(result$dissimilarity - result_base_cor)),
    1e-9
  )
})


# =============================================================================
# diss_cosine dissimilarity tests
# =============================================================================

test_that("diss_cosine works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  result <- dissimilarity(
    Xr = d$Xr, 
    Xu = d$Xu,
    diss_method = diss_cosine(center = TRUE, scale = FALSE)
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
  expect_equal(nrow(result$dissimilarity), nrow(d$Xr))
  expect_equal(ncol(result$dissimilarity), nrow(d$Xu))
})


test_that("diss_cosine self-dissimilarity diagonal is zero", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 20, n_xu = 0)
  
  result <- dissimilarity(
    Xr = d$Xr,
    diss_method = diss_cosine()
  )
  
  expect_true(all(abs(diag(result$dissimilarity)) < 1e-6))
})


# =============================================================================
# diss_mahalanobis dissimilarity tests
# =============================================================================

test_that("diss_mahalanobis works when n > p", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 40, n_xu = 20)
  
  # Use only first 20 columns to ensure n > p
  result <- dissimilarity(
    Xr = d$Xr[, 1:20], 
    Xu = d$Xu[, 1:20],
    diss_method = diss_mahalanobis(center = TRUE, scale = FALSE)
  )
  
  expect_s3_class(result, "dissimilarity")
  expect_true("dissimilarity" %in% names(result))
})


test_that("diss_mahalanobis errors when covariance is singular", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 10, n_xu = 5)
  
  # Full spectra: p >> n, covariance will be singular
  expect_error(
    dissimilarity(
      Xr = d$Xr,
      Xu = d$Xu,
      diss_method = diss_mahalanobis()
    ),
    "number of observations.*larger than the number of variables"
  )
})


# =============================================================================
# Legacy API tests
# =============================================================================

test_that("character-based diss_method is rejected with helpful message", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data()
  
  expect_error(
    dissimilarity(Xr = d$Xr, Xu = d$Xu, diss_method = "pca"),
    "no longer supported"
  )
  
  expect_error(
    dissimilarity(Xr = d$Xr, Xu = d$Xu, diss_method = "euclid"),
    "diss_euclidean"
  )
})


# =============================================================================
# Large dataset tests (skipped on CRAN)
# =============================================================================

test_that("dissimilarity with large datasets works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  Xr_snv <- prospectr::standardNormalVariate(Xr)
  Xu_snv <- prospectr::standardNormalVariate(Xu)
  
  # PCA with opc
  result_pca <- dissimilarity(
    Xr = Xr, 
    Xu = Xu,
    Yr = Yr,
    diss_method = diss_pca(
      ncomp = ncomp_by_opc(max_ncomp = 30),
      center = TRUE, 
      scale = TRUE,
      return_projection = TRUE
    )
  )
  
  # PLS with opc
  result_pls <- dissimilarity(
    Xr = Xr_snv, 
    Xu = Xu_snv,
    Yr = Yr,
    diss_method = diss_pls(
      ncomp = ncomp_by_opc(max_ncomp = 30),
      return_projection = TRUE
    )
  )
  
  # Euclidean
  result_euclid <- dissimilarity(
    Xr = Xr, 
    Xu = Xu,
    diss_method = diss_euclidean(center = TRUE, scale = TRUE)
  )
  
  # Correlation with moving window
  result_cor <- dissimilarity(
    Xr = Xr, 
    Xu = Xu,
    diss_method = diss_correlation(ws = 11)
  )
  
  expect_s3_class(result_pca, "dissimilarity")
  expect_s3_class(result_pls, "dissimilarity")
  expect_s3_class(result_euclid, "dissimilarity")
  expect_s3_class(result_cor, "dissimilarity")
  
  expect_true(result_pca$ncomp <= 30)
  expect_true(result_pls$ncomp <= 30)
  
  # Verify Euclidean against base R
  result_euclid_self <- dissimilarity(
    Xr = Xu[1:10, ],
    diss_method = diss_euclidean(center = FALSE, scale = FALSE)
  )
  resemble_diss <- result_euclid_self$dissimilarity
  resemble_diss_unscaled <- ((resemble_diss^2) * ncol(Xu))^0.5
  base_diss <- as.matrix(dist(Xu[1:10, ]))
  
  expect_lt(sum(abs(resemble_diss_unscaled - base_diss)), 1e-5)
})


# =============================================================================
# diss_evaluate tests
# =============================================================================

test_that("diss_evaluate works with continuous side_info", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 40, n_xu = 0)
  
  result <- dissimilarity(
    Xr = d$Xr,
    diss_method = diss_pca(ncomp = 8)
  )
  
  ev <- diss_evaluate(
    diss = result$dissimilarity, 
    side_info = as.matrix(d$Yr)
  )
  
  expect_true("eval" %in% names(ev))
  expect_true("first_nn" %in% names(ev))
  expect_true("rmsd" %in% colnames(ev$eval))
  expect_true("r" %in% colnames(ev$eval))
})


test_that("diss_evaluate works with multiple side_info columns", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Yr2 <- NIRsoil$Ciso[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  complete <- !is.na(Yr) & !is.na(Yr2)
  Xr <- Xr[complete, ][1:40, ]
  Yr <- Yr[complete][1:40]
  Yr2 <- Yr2[complete][1:40]
  
  result <- dissimilarity(
    Xr = Xr,
    diss_method = diss_pca(ncomp = 8)
  )
  
  ev <- diss_evaluate(
    diss = result$dissimilarity, 
    side_info = cbind(Yr, Yr2)
  )
  
  expect_true("eval" %in% names(ev))
  expect_true("global_eval" %in% names(ev))
  expect_equal(nrow(ev$eval), 2)
})


test_that("diss_evaluate validates inputs", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_nirsoil_data(n_xr = 20, n_xu = 0)
  
  result <- dissimilarity(
    Xr = d$Xr,
    diss_method = diss_pca(ncomp = 5)
  )
  
  # side_info must be a matrix
  expect_error(
    diss_evaluate(result$dissimilarity, side_info = d$Yr),
    "must be a matrix"
  )
  
  # Dimension mismatch
  expect_error(
    diss_evaluate(result$dissimilarity, side_info = as.matrix(d$Yr[1:10])),
    "must match"
  )
})





# =============================================================================
# Tests for cor_diss (internal correlation dissimilarity function)
# =============================================================================

# -----------------------------------------------------------------------------
# Basic functionality tests
# -----------------------------------------------------------------------------

test_that("cor_diss returns correct dimensions for Xr only", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- cor_diss(d$Xr)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xr))
})

test_that("cor_diss returns correct dimensions with Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- cor_diss(d$Xr, d$Xu)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xu))
})

test_that("cor_diss diagonal is zero for self-dissimilarity", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- cor_diss(d$Xr)
  
  expect_true(all(abs(diag(result)) < 1e-6))
})

test_that("cor_diss is symmetric", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- cor_diss(d$Xr)
  
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("cor_diss values are bounded [0, 1]", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- cor_diss(d$Xr)
  
  expect_true(all(result >= -1e-10))
  expect_true(all(result <= 1 + 1e-10))
})

test_that("cor_diss with moving window returns correct dimensions", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- cor_diss(d$Xr, ws = 41)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xr))
})

test_that("cor_diss with moving window and Xu works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- cor_diss(d$Xr, d$Xu, ws = 41)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xu))
})

test_that("cor_diss with different ws values gives different results", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result_no_ws <- cor_diss(d$Xr)
  result_ws41 <- cor_diss(d$Xr, ws = 41)
  result_ws21 <- cor_diss(d$Xr, ws = 21)
  
  expect_false(all(result_no_ws == result_ws41))
  expect_false(all(result_ws41 == result_ws21))
})

test_that("cor_diss centering and scaling work", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result_centered <- cor_diss(d$Xr, center = TRUE)
  result_scaled <- cor_diss(d$Xr, scale = TRUE)
  result_both <- cor_diss(d$Xr, center = TRUE, scale = TRUE)
  
  expect_true(is.matrix(result_centered))
  expect_true(is.matrix(result_scaled))
  expect_true(is.matrix(result_both))
})

test_that("cor_diss row/column names are correct", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data(n_xr = 5, n_xu = 3)
  
  result_xr <- cor_diss(d$Xr)
  result_xu <- cor_diss(d$Xr, d$Xu)
  
  expect_equal(rownames(result_xr), paste0("Xr_", 1:5))
  expect_equal(colnames(result_xr), paste0("Xr_", 1:5))
  expect_equal(rownames(result_xu), paste0("Xr_", 1:5))
  expect_equal(colnames(result_xu), paste0("Xu_", 1:3))
})

# -----------------------------------------------------------------------------
# Error handling tests
# -----------------------------------------------------------------------------

test_that("cor_diss errors with single variable", {
  skip_on_cran()
  Xr <- matrix(1:10, ncol = 1)
  
  expect_error(cor_diss(Xr), "number of variables must be larger than 1")
})

test_that("cor_diss errors with mismatched columns", {
  skip_on_cran()
  Xr <- matrix(1:20, nrow = 2)
  Xu <- matrix(1:10, nrow = 2)
  
  expect_error(cor_diss(Xr, Xu), "number of columns")
})

test_that("cor_diss errors with NA values", {
  skip_on_cran()
  Xr <- matrix(1:20, nrow = 2)
  Xr[1, 1] <- NA
  
  expect_error(cor_diss(Xr), "missing values")
  
  Xr <- matrix(1:20, nrow = 2)
  Xu <- matrix(1:20, nrow = 2)
  Xu[1, 1] <- NA
  
  expect_error(cor_diss(Xr, Xu), "missing values")
})

test_that("cor_diss errors with invalid logical arguments", {
  skip_on_cran()
  Xr <- matrix(1:20, nrow = 2)
  
  expect_error(cor_diss(Xr, center = "yes"), "'center' argument must be logical")
  expect_error(cor_diss(Xr, scale = 1), "'scale' argument must be logical")
})

test_that("cor_diss errors with invalid ws", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  expect_error(cor_diss(d$Xr, ws = 40), "odd value")
  expect_error(cor_diss(d$Xr, ws = 1), "greater than 2")
  expect_error(cor_diss(d$Xr, ws = ncol(d$Xr) + 1), "smaller than the number of columns")
})

test_that("cor_diss errors with zero standard deviation observation", {
  skip_on_cran()
  Xr <- matrix(rnorm(100), nrow = 10)
  Xr[1, ] <- rep(5, 10)
  
  expect_error(cor_diss(Xr, center = FALSE), "Correlation coefficients cannot be computed")
})

# =============================================================================
# Tests for .f_diss_compute (internal function)
# =============================================================================

test_that(".f_diss_compute works with euclidean method", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  method <- diss_euclidean(center = TRUE, scale = FALSE)
  result <- .f_diss_compute(d$Xr, Xu = NULL, method = method)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xr))
  expect_true(all(abs(diag(result)) < 1e-6))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that(".f_diss_compute works with euclidean and Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  method <- diss_euclidean()
  result <- .f_diss_compute(d$Xr, Xu = d$Xu, method = method)
  
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xu))
  expect_equal(rownames(result), paste0("Xr_", seq_len(nrow(d$Xr))))
  expect_equal(colnames(result), paste0("Xu_", seq_len(nrow(d$Xu))))
})

test_that(".f_diss_compute works with euclidean center and scale", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  method_none <- diss_euclidean(center = FALSE, scale = FALSE)
  method_center <- diss_euclidean(center = TRUE, scale = FALSE)
  method_scale <- diss_euclidean(center = FALSE, scale = TRUE)
  method_both <- diss_euclidean(center = TRUE, scale = TRUE)
  
  result_none <- .f_diss_compute(d$Xr, method = method_none)
  result_center <- .f_diss_compute(d$Xr, method = method_center)
  result_scale <- .f_diss_compute(d$Xr, method = method_scale)
  result_both <- .f_diss_compute(d$Xr, method = method_both)
  
  # All should produce valid matrices
  expect_true(is.matrix(result_none))
  expect_true(is.matrix(result_center))
  expect_true(is.matrix(result_scale))
  expect_true(is.matrix(result_both))
  
  # Results should differ
  expect_false(all(result_none == result_center))
  expect_false(all(result_none == result_scale))
})

test_that(".f_diss_compute works with mahalanobis method", {
  skip_on_cran()
  
  # Create data where n > p (required for Mahalanobis)
  set.seed(123)
  Xr <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  
  method <- diss_mahalanobis(center = TRUE, scale = FALSE)
  result <- .f_diss_compute(Xr, Xu = NULL, method = method)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(Xr))
  expect_equal(ncol(result), nrow(Xr))
  expect_true(all(abs(diag(result)) < 1e-6))
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that(".f_diss_compute works with mahalanobis and Xu", {
  skip_on_cran()
  
  set.seed(123)
  Xr <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  Xu <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
  
  method <- diss_mahalanobis()
  result <- .f_diss_compute(Xr, Xu = Xu, method = method)
  
  expect_equal(nrow(result), nrow(Xr))
  expect_equal(ncol(result), nrow(Xu))
})

test_that(".f_diss_compute errors with singular covariance for mahalanobis", {
  skip_on_cran()
  
  # Create data where n < p (singular covariance)
  set.seed(123)
  Xr <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
  
  method <- diss_mahalanobis()
  
  expect_error(
    .f_diss_compute(Xr, method = method),
    "singular"
  )
})

test_that(".f_diss_compute works with cosine method", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  method <- diss_cosine(center = FALSE, scale = FALSE)
  result <- .f_diss_compute(d$Xr, Xu = NULL, method = method)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xr))
  expect_true(all(abs(diag(result)) < 1e-6))
  expect_equal(result, t(result), tolerance = 1e-10)
  
  # Cosine dissimilarity should be bounded [0, pi]
  expect_true(all(result >= -1e-10))
  expect_true(all(result <= pi + 1e-10))
})

test_that(".f_diss_compute works with cosine and Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  method <- diss_cosine()
  result <- .f_diss_compute(d$Xr, Xu = d$Xu, method = method)
  
  expect_equal(nrow(result), nrow(d$Xr))
  expect_equal(ncol(result), nrow(d$Xu))
})

test_that(".f_diss_compute cosine handles identical observations", {
  skip_on_cran()
  
  # Two identical rows should have zero dissimilarity
  Xr <- matrix(c(1, 2, 3, 1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
  
  method <- diss_cosine(center = FALSE, scale = FALSE)
  result <- .f_diss_compute(Xr, method = method)
  
  # First two rows are identical - should be zero
  expect_equal(result[1, 2], 0, tolerance = 1e-10)
  expect_equal(result[2, 1], 0, tolerance = 1e-10)
})

test_that(".f_diss_compute cosine replaces NaN with zero", {
  skip_on_cran()
  
  # Perfectly identical rows can cause NaN in acos due to floating point
  Xr <- matrix(rep(1:10, 3), nrow = 3, byrow = TRUE)
  
  method <- diss_cosine(center = FALSE, scale = FALSE)
  result <- .f_diss_compute(Xr, method = method)
  
  # Should have no NaN values
  expect_false(any(is.nan(result)))
  
  # Diagonal and identical rows should be zero
  expect_true(all(abs(diag(result)) < 1e-10))
})

test_that(".f_diss_compute cosine with center and scale", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  method_none <- diss_cosine(center = FALSE, scale = FALSE)
  method_center <- diss_cosine(center = TRUE, scale = FALSE)
  method_both <- diss_cosine(center = TRUE, scale = TRUE)
  
  result_none <- .f_diss_compute(d$Xr, method = method_none)
  result_center <- .f_diss_compute(d$Xr, method = method_center)
  result_both <- .f_diss_compute(d$Xr, method = method_both)
  
  expect_true(is.matrix(result_none))
  expect_true(is.matrix(result_center))
  expect_true(is.matrix(result_both))
})

test_that(".f_diss_compute row/column names are correct", {
  skip_on_cran()
  
  Xr <- matrix(rnorm(25), nrow = 5)
  Xu <- matrix(rnorm(15), nrow = 3)
  
  method <- diss_euclidean()
  
  result_xr <- .f_diss_compute(Xr, method = method)
  expect_equal(rownames(result_xr), paste0("Xr_", 1:5))
  expect_equal(colnames(result_xr), paste0("Xr_", 1:5))
  
  result_xu <- .f_diss_compute(Xr, Xu = Xu, method = method)
  expect_equal(rownames(result_xu), paste0("Xr_", 1:5))
  expect_equal(colnames(result_xu), paste0("Xu_", 1:3))
})

test_that(".f_diss_compute validates method class", {
  skip_on_cran()
  
  Xr <- matrix(rnorm(20), nrow = 4)
  
  expect_error(.f_diss_compute(Xr, method = "euclidean"))
  expect_error(.f_diss_compute(Xr, method = list(method = "euclidean")))
})

test_that(".f_diss_compute validates method name", {
  skip_on_cran()
  
  Xr <- matrix(rnorm(20), nrow = 4)
  
  # Create a fake diss_method with invalid method name
  fake_method <- structure(
    list(method = "invalid", center = TRUE, scale = FALSE),
    class = c("diss_invalid", "diss_method")
  )
  
  expect_error(.f_diss_compute(Xr, method = fake_method))
})