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
# ncomp_* constructor tests
# =============================================================================

test_that("ncomp_by_var constructor works", {
  nc <- ncomp_by_var(0.01)
  expect_s3_class(nc, c("ncomp_by_var", "ncomp_selection"))
  expect_equal(nc$min_var, 0.01)
  expect_equal(nc$max_ncomp, 40L)
  
  nc2 <- ncomp_by_var(0.05, max_ncomp = 20)
  expect_equal(nc2$min_var, 0.05)
  expect_equal(nc2$max_ncomp, 20L)
})


test_that("ncomp_by_var validates inputs", {
  expect_error(ncomp_by_var(0), "\\(0, 1\\]")
  expect_error(ncomp_by_var(1.5), "\\(0, 1\\]")
  expect_error(ncomp_by_var("a"), "\\(0, 1\\]")
  expect_error(ncomp_by_var(0.01, max_ncomp = 0), "positive integer")
})


test_that("ncomp_by_cumvar constructor works", {
  nc <- ncomp_by_cumvar(0.99)
  expect_s3_class(nc, c("ncomp_by_cumvar", "ncomp_selection"))
  expect_equal(nc$min_cumvar, 0.99)
  expect_equal(nc$max_ncomp, 40L)
  
  nc2 <- ncomp_by_cumvar(0.95, max_ncomp = 50)
  expect_equal(nc2$min_cumvar, 0.95)
  expect_equal(nc2$max_ncomp, 50L)
})


test_that("ncomp_by_cumvar validates inputs", {
  expect_error(ncomp_by_cumvar(0), "\\(0, 1\\]")
  expect_error(ncomp_by_cumvar(1.5), "\\(0, 1\\]")
  expect_error(ncomp_by_cumvar(0.99, max_ncomp = -1), "positive integer")
})


test_that("ncomp_by_opc constructor works", {
  nc <- ncomp_by_opc()
  expect_s3_class(nc, c("ncomp_by_opc", "ncomp_selection"))
  expect_equal(nc$max_ncomp, 40L)
  
  nc2 <- ncomp_by_opc(max_ncomp = 30)
  expect_equal(nc2$max_ncomp, 30L)
})


test_that("ncomp_by_opc validates inputs", {
  expect_error(ncomp_by_opc(max_ncomp = 0), "positive integer")
  expect_error(ncomp_by_opc(max_ncomp = "a"), "positive integer")
})


test_that("ncomp_fixed constructor works", {
  nc <- ncomp_fixed(10)
  expect_s3_class(nc, c("ncomp_fixed", "ncomp_selection"))
  expect_equal(nc$ncomp, 10L)
})


test_that("ncomp_fixed validates inputs", {
  expect_error(ncomp_fixed(), "required")
  expect_error(ncomp_fixed(0), "positive integer")
  expect_error(ncomp_fixed(-5), "positive integer")
  expect_error(ncomp_fixed("a"), "positive integer")
})


# =============================================================================
# diss_* constructor tests
# =============================================================================

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
