context("test-search_neighbors")
# tests/testthat/test-search_neighbors.R
library(foreach)
library(RhpcBLASctl)
registerDoSEQ()
# =============================================================================
# Setup helper
# =============================================================================

.setup_search_neighbors_data <- function() {
  data("NIRsoil", package = "prospectr")
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  list(Xr = Xr, Xu = Xu, Yr = Yr, Yu = Yu)
}


# =============================================================================
# Constructor validation tests
# =============================================================================

test_that("neighbors_k constructor works", {
  n1 <- neighbors_k(40)
  expect_s3_class(n1, "neighbors_k")
  expect_s3_class(n1, "neighbors")
  expect_equal(n1$k, 40L)
  
  # Multiple k values
  n2 <- neighbors_k(c(10, 20, 30))
  expect_equal(n2$k, c(10L, 20L, 30L))
  
  # Validation: k must be positive
  expect_error(neighbors_k(0))
  expect_error(neighbors_k(-5))
})


test_that("neighbors_diss constructor works", {
  n1 <- neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  expect_s3_class(n1, "neighbors_diss")
  expect_s3_class(n1, "neighbors")
  expect_equal(n1$threshold, 0.1)
  expect_equal(n1$k_min, 10L)
  expect_equal(n1$k_max, 50L)
  
  # Default k_max is Inf
  n2 <- neighbors_diss(threshold = 0.5)
  expect_equal(n2$k_max, Inf)
})


# =============================================================================
# Legacy API rejection tests
# =============================================================================

test_that("search_neighbors rejects legacy arguments", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  # Reject old k argument
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = diss_pca(),
      k = 50,
      neighbors = neighbors_k(50)
    ),
    "k.*k_diss.*k_range.*removed"
  )
  
  # Reject old pc_selection argument
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = diss_pca(),
      neighbors = neighbors_k(50),
      pc_selection = list("manual", 10)
    ),
    "pc_selection.*removed"
  )
  
  # Reject old center/scale arguments
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = diss_pca(),
      neighbors = neighbors_k(50),
      center = TRUE
    ),
    "center.*scale.*removed"
  )
  
  # Reject character-based diss_method
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = "pca",
      neighbors = neighbors_k(50)
    ),
    "Character-based.*no longer supported"
  )
})


test_that("search_neighbors requires neighbors argument", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = diss_pca()
    ),
    "neighbors.*required"
  )
})


# =============================================================================
# Basic functionality tests
# =============================================================================

test_that("search_neighbors with diss_pca and neighbors_k works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10, method = "pca_nipals"),
    neighbors = neighbors_k(50)
  )
  
  expect_type(nn, "list")
  expect_true("neighbors" %in% names(nn))
  expect_true("neighbors_diss" %in% names(nn))
  expect_true("unique_neighbors" %in% names(nn))
  expect_equal(nrow(nn$neighbors), 50)
  expect_equal(ncol(nn$neighbors), nrow(d$Xu))
})


test_that("search_neighbors with diss_pca and neighbors_diss works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10, method = "pca_nipals"),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  )
  
  expect_type(nn, "list")
  expect_true("neighbors" %in% names(nn))
  expect_true("k_diss_info" %in% names(nn))
})


test_that("search_neighbors with spike works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  spikes <- c(5, 10, 15)
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10, method = "pca_nipals"),
    neighbors = neighbors_k(50),
    spike = spikes
  )
  
  # Spiked observations should appear in all neighborhoods
  expect_true(all(rowMeans(nn$neighbors[seq_along(spikes), ]) == spikes))
})


test_that("search_neighbors with spike and neighbors_diss works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  spikes <- c(5, 10, 15)
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10, method = "pca_nipals"),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50),
    spike = spikes
  )
  
  # Spiked observations should appear in all neighborhoods
  expect_true(all(rowMeans(nn$neighbors[seq_along(spikes), ]) == spikes))
})


# =============================================================================
# Dissimilarity method tests
# =============================================================================

test_that("search_neighbors works with diss_pca (SVD)", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10),
    Yr = d$Yr,
    neighbors = neighbors_k(10)
  )
  
  expect_type(nn, "list")
  expect_equal(nrow(nn$neighbors), 10)
})


test_that("search_neighbors works with diss_pls", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pls(ncomp = 10),
    Yr = d$Yr,
    neighbors = neighbors_k(10)
  )
  
  expect_type(nn, "list")
  expect_equal(nrow(nn$neighbors), 10)
})


test_that("search_neighbors works with diss_correlation", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_correlation(),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  )
  
  expect_type(nn, "list")
})


test_that("search_neighbors works with diss_euclidean", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_euclidean(),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  )
  
  expect_type(nn, "list")
})


test_that("search_neighbors works with diss_cosine", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_cosine(),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  )
  
  expect_type(nn, "list")
})


# =============================================================================
# Return options tests
# =============================================================================

test_that("search_neighbors with return_dissimilarity works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(50),
    return_dissimilarity = TRUE
  )
  
  expect_true("dissimilarity" %in% names(nn))
})


test_that("search_neighbors with return_projection works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(
      ncomp = ncomp_by_opc(40),
      return_projection = TRUE
    ),
    Yr = d$Yr,
    neighbors = neighbors_k(50)
  )
  
  expect_true("projection" %in% names(nn))
})


# =============================================================================
# Expected results tests (skip on CRAN)
# =============================================================================

test_that("search_neighbors produces expected neighbor indices", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  # PCA-based
  nn_pca <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(10)
  )
  expected_sum_pca <- 17853
  expect_equal(sum(nn_pca$neighbors[1, ]), expected_sum_pca)
  
  # PLS-based
  nn_pls <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_pls(ncomp = 10),
    Yr = d$Yr,
    neighbors = neighbors_k(10)
  )
  expected_sum_pls <- 18410
  expect_equal(sum(nn_pls$neighbors[1, ]), expected_sum_pls)
  
  # Correlation-based
  nn_cor <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_correlation(),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  )
  expected_sum_cor <- 18036
  expect_equal(sum(nn_cor$neighbors[1, ]), expected_sum_cor)
  
  # Euclidean
  nn_euclid <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_euclidean(),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  )
  expected_sum_euclid <- 19586
  expect_equal(sum(nn_euclid$neighbors[1, ]), expected_sum_euclid)
  
  # Cosine
  nn_cosine <- search_neighbors(
    Xr = d$Xr, Xu = d$Xu,
    diss_method = diss_cosine(),
    neighbors = neighbors_diss(threshold = 0.1, k_min = 10, k_max = 50)
  )
  expected_sum_cosine <- 18910
  expect_equal(sum(nn_cosine$neighbors[1, ]), expected_sum_cosine)
})


# =============================================================================
# Edge cases and validation tests
# =============================================================================

test_that("search_neighbors validates k vs nrow(Xr)", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = diss_pca(ncomp = 10),
      neighbors = neighbors_k(nrow(d$Xr) + 10)
    ),
    "cannot exceed"
  )
})


test_that("search_neighbors validates spike indices", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  # Out of bounds
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = diss_pca(ncomp = 10),
      neighbors = neighbors_k(50),
      spike = c(1, nrow(d$Xr) + 100)
    ),
    "out-of-bounds"
  )
  
  # Contradictory indices
  expect_error(
    search_neighbors(
      Xr = d$Xr, Xu = d$Xu,
      diss_method = diss_pca(ncomp = 10),
      neighbors = neighbors_k(50),
      spike = c(5, -5)
    ),
    "contradictory"
  )
})


test_that("search_neighbors without Xu searches within Xr", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_search_neighbors_data()
  
  nn <- search_neighbors(
    Xr = d$Xr,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(20)
  )
  
  expect_type(nn, "list")
  expect_equal(ncol(nn$neighbors), nrow(d$Xr))
  
  # Column names should be Xr* not Xu*
  expect_true(all(grepl("^Xr", colnames(nn$neighbors))))
})

