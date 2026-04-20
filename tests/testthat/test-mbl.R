context("test-mbl")
# tests/testthat/test-mbl.R
library(foreach)
library(RhpcBLASctl)
registerDoSEQ()
# =============================================================================
# Setup helper
# =============================================================================

.setup_mbl_data <- function(n_xr = 40, n_xu = 20, preprocess = FALSE) {
  data("NIRsoil", package = "prospectr")
  
  if (preprocess) {
    NIRsoil$spc <- prospectr::savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
  }
  
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
# mbl_control tests
# =============================================================================

test_that("mbl_control constructor works", {
  ctrl <- mbl_control()
  
  expect_type(ctrl, "list")
  expect_false(ctrl$return_dissimilarity)
  expect_equal(ctrl$validation_type, "NNv")
  expect_true(ctrl$tune_locally)
  expect_equal(ctrl$number, 10L)
  expect_equal(ctrl$p, 0.75)
  expect_true(ctrl$range_prediction_limits)
  expect_true(ctrl$allow_parallel)
  expect_equal(ctrl$blas_threads, 1L)
})


test_that("mbl_control validates inputs", {
  expect_error(mbl_control(return_dissimilarity = "yes"), "TRUE or FALSE")
  expect_error(
    mbl_control(validation_type = "invalid"),
    "NNv.*local_cv.*none"
  )
  expect_error(mbl_control(tune_locally = "yes"), "TRUE or FALSE")
  expect_error(mbl_control(validation_type = "local_cv", number = 0), "positive integer")
  expect_error(mbl_control(validation_type = "local_cv", p = 1.5), "between 0 and 1")
  expect_error(mbl_control(validation_type = "local_cv", p = 0), "between 0 and 1")
  expect_error(mbl_control(allow_parallel = "yes"), "TRUE or FALSE")
  expect_error(mbl_control(blas_threads = 0), "positive integer")
})


test_that("mbl_control rejects 'none' combined with other types", {
  expect_error(
    mbl_control(validation_type = c("NNv", "none")),
    "cannot combine 'none'"
  )
})


# =============================================================================
# Basic mbl functionality tests
# =============================================================================

test_that("mbl with fit_gpr and neighbors_k works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.5,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_gpr(),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
  expect_true("validation_results" %in% names(result))
  expect_true("results" %in% names(result))
  expect_true("Xu_neighbors" %in% names(result))
})


test_that("mbl with fit_pls and neighbors_k works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.5,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
  expect_true("validation_results" %in% names(result))
})


test_that("mbl with fit_pls (mpls method) works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_pls(ncomp = 5, method = "mpls"),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with fit_pls (simpls method) works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_pls(ncomp = 5, method = "simpls"),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with fit_wapls and neighbors_k works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.5,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with fit_wapls (mpls method) works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 5, method = "mpls"),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with fit_wapls (simpls method) works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 5, method = "simpls"),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


# =============================================================================
# neighbors_diss tests
# =============================================================================

test_that("mbl with fit_gpr and neighbors_diss works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.5,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_diss(threshold = 0.1, k_min = 15, k_max = 30),
    fit_method = fit_gpr(),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with fit_pls and neighbors_diss works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_diss(threshold = 0.1, k_min = 15, k_max = 30),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with fit_wapls and neighbors_diss works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_diss(threshold = 0.1, k_min = 15, k_max = 30),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


# =============================================================================
# Group argument tests
# =============================================================================

test_that("mbl with group argument works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.5,
    allow_parallel = FALSE
  )
  
  group <- rep(c(1, 2), length.out = nrow(d$Xr))
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    group = group,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


# =============================================================================
# Dissimilarity method tests
# =============================================================================

test_that("mbl with diss_pca works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_pca(ncomp = 10),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with diss_pls works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_pls(ncomp = 10),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with diss_correlation works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_correlation(ws = 11),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with diss_euclidean works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_euclidean(center = TRUE, scale = TRUE),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with diss_cosine works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_cosine(),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


# =============================================================================
# External dissimilarity matrix tests
# =============================================================================

test_that("mbl with external dissimilarity matrix works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  tol <- 1e-8
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  # Compute external dissimilarity matrix
  combined <- rbind(d$Xr, d$Xu)
  ext_diss <- dissimilarity(
    Xr = combined,
    diss_method = diss_correlation(center = FALSE, scale = FALSE)
  )$dissimilarity
  diag(ext_diss) <- 0
  
  # mbl with external matrix
  result_ext <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = ext_diss,
    diss_usage = "predictors",
    fit_method = fit_gpr(),
    control = ctrl,
    verbose = FALSE
  )
  
  # mbl with internal computation (same method)
  result_int <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_correlation(center = FALSE, scale = FALSE),
    diss_usage = "predictors",
    fit_method = fit_gpr(),
    control = ctrl,
    verbose = FALSE
  )
  
  # Results should be identical
  nnv_ext <- result_ext$validation_results$nearest_neighbor_validation
  nnv_int <- result_int$validation_results$nearest_neighbor_validation
  
  expect_lt(sum(abs(nnv_ext - nnv_int)), tol)
})


# =============================================================================
# diss_usage tests
# =============================================================================

test_that("mbl with diss_usage = 'predictors' works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_pca(ncomp = 10),
    diss_usage = "predictors",
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


test_that("mbl with diss_usage = 'weights' works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    diss_method = diss_pca(ncomp = 10),
    diss_usage = "weights",
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
})


# =============================================================================
# spike argument tests
# =============================================================================

test_that("mbl with spike argument works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_gpr(),
    spike = 1:5,
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
  expect_true(result$spike)
})


# =============================================================================
# gh argument tests
# =============================================================================

test_that("mbl with gh = TRUE works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(c(25, 35)),
    fit_method = fit_pls(ncomp = 5),
    gh = TRUE,
    control = ctrl,
    verbose = FALSE
  )
  
  expect_s3_class(result, "mbl")
  expect_true(!is.null(result$gh))
})


# =============================================================================
# Expected results tests (skipped on CRAN)
# =============================================================================
## it's a sanity check ensuring results stay within plausible bounds rather 
## than testing for exact values. This catches regressions where something 
## breaks catastrophically 
test_that("mbl delivers expected results", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  NIRsoil$spc <- prospectr::savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.8,
    tune_locally = TRUE,
    allow_parallel = FALSE
  )
  
  k_test <- c(40, 150)
  tseed <- 141020
  
  # GPR
  set.seed(tseed)
  gpr <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_k(k_test),
    fit_method = fit_gpr(noise_variance = 0.0001, scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  # PLS
  set.seed(tseed)
  pls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_k(k_test),
    fit_method = fit_pls(ncomp = 10, scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  # WAPLS
  set.seed(tseed)
  wapls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_k(k_test),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15, scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  # MPLS
  set.seed(tseed)
  mpls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_k(k_test),
    fit_method = fit_pls(ncomp = 10, method = "mpls", scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  # WAMPLS
  set.seed(tseed)
  wampls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_k(k_test),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15, method = "mpls", scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  
  wampls_simpls <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_k(k_test),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15, method = "simpls", scale = TRUE),
    control = ctrl,
    verbose = FALSE
  )
  
  # Check local CV RMSE bounds
  expect_true(all(
    gpr$validation_results$local_cross_validation$rmse < 2.5 &
      gpr$validation_results$local_cross_validation$rmse > 1.5
  ))
  
  expect_true(all(
    pls$validation_results$local_cross_validation$rmse < 2 &
      pls$validation_results$local_cross_validation$rmse > 1.4
  ))
  
  expect_true(all(
    wapls$validation_results$local_cross_validation$rmse < 1.8 &
      wapls$validation_results$local_cross_validation$rmse > 1.5
  ))
  
  expect_true(all(
    mpls$validation_results$local_cross_validation$rmse < 2 &
      mpls$validation_results$local_cross_validation$rmse > 1.5
  ))
  
  expect_true(all(
    wampls$validation_results$local_cross_validation$rmse < 1.9 &
      wampls$validation_results$local_cross_validation$rmse > 1.5
  ))
  # WAMPLS SIMPLS bounds (should be similar to WAPLS/WAMPLS)
  expect_true(all(
    wampls_simpls$validation_results$local_cross_validation$rmse < 1.9 &
      wampls_simpls$validation_results$local_cross_validation$rmse > 1.5
  ))
  
  
  # Check NNv R2 bounds
  expect_true(all(gpr$validation_results$nearest_neighbor_validation$r2 > 0.50))
  expect_true(all(pls$validation_results$nearest_neighbor_validation$r2 > 0.74))
  expect_true(all(wapls$validation_results$nearest_neighbor_validation$r2 > 0.80))
  
  # Check Yu prediction R2 bounds
  expect_true(all(gpr$validation_results$Yu_prediction_statistics$r2 > 0.67))
  expect_true(all(pls$validation_results$Yu_prediction_statistics$r2 > 0.60))
  expect_true(all(wapls$validation_results$Yu_prediction_statistics$r2 > 0.69))
  expect_true(all(wampls_simpls$validation_results$nearest_neighbor_validation$r2 > 0.78))
  expect_true(all(wampls_simpls$validation_results$Yu_prediction_statistics$r2 > 0.65))
})


test_that("mbl with neighbors_diss delivers expected results", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  NIRsoil$spc <- prospectr::savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.8,
    tune_locally = TRUE,
    allow_parallel = FALSE
  )
  
  tseed <- 141020
  
  # GPR with neighbors_diss
  set.seed(tseed)
  gpr_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_diss(threshold = 0.1, k_min = 20, k_max = 100),
    fit_method = fit_gpr(),
    control = ctrl,
    verbose = FALSE
  )
  
  # PLS with neighbors_diss
  set.seed(tseed)
  pls_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_diss(threshold = 0.1, k_min = 20, k_max = 100),
    fit_method = fit_pls(ncomp = 10),
    control = ctrl,
    verbose = FALSE
  )
  
  # WAPLS with neighbors_diss
  set.seed(tseed)
  wapls_diss <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_diss(threshold = 0.1, k_min = 20, k_max = 100),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 15),
    control = ctrl,
    verbose = FALSE
  )
  
  # Check bounds
  expect_true(all(
    gpr_diss$validation_results$local_cross_validation$rmse < 2.8 &
      gpr_diss$validation_results$local_cross_validation$rmse > 2.5
  ))
  
  expect_true(all(
    pls_diss$validation_results$local_cross_validation$rmse < 1.8 &
      pls_diss$validation_results$local_cross_validation$rmse > 1.4
  ))
  
  expect_true(all(
    wapls_diss$validation_results$local_cross_validation$rmse < 2 &
      wapls_diss$validation_results$local_cross_validation$rmse > 1.4
  ))
  
  expect_true(all(gpr_diss$validation_results$nearest_neighbor_validation$r2 > 0.76))
  expect_true(all(pls_diss$validation_results$nearest_neighbor_validation$r2 > 0.81))
  expect_true(all(wapls_diss$validation_results$nearest_neighbor_validation$r2 > 0.81))
  
  expect_true(all(gpr_diss$validation_results$Yu_prediction_statistics$r2 > 0.60))
  expect_true(all(pls_diss$validation_results$Yu_prediction_statistics$r2 > 0.60))
  expect_true(all(wapls_diss$validation_results$Yu_prediction_statistics$r2 > 0.65))
})


test_that("mbl with group argument delivers expected results", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  NIRsoil$spc <- prospectr::savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0)
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4, 
    p = 0.8,
    tune_locally = TRUE,
    allow_parallel = FALSE
  )
  
  k_test <- c(40, 150)
  tseed <- 141020
  
  xgroup <- rep(seq_len(floor(nrow(Xr) / 2)), each = 2)
  
  set.seed(tseed)
  pls_group <- mbl(
    Xr = Xr, Yr = Yr, Xu = Xu, Yu = Yu,
    neighbors = neighbors_k(k_test),
    fit_method = fit_pls(ncomp = 10),
    control = ctrl, 
    group = xgroup,
    verbose = FALSE
  )
  
  expect_true(all(
    pls_group$validation_results$local_cross_validation$rmse < 2 &
      pls_group$validation_results$local_cross_validation$rmse > 1.4
  ))
  
  expect_true(all(pls_group$validation_results$nearest_neighbor_validation$r2 > 0.7))
})


.with_test_device <- function(code) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)
  force(code)
}

.build_mbl_for_plot <- function(gh = FALSE) {
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    allow_parallel = FALSE
  )
  
  mbl(
    Xr = d$Xr,
    Yr = d$Yr,
    Xu = d$Xu,
    Yu = d$Yu,
    neighbors = neighbors_k(25),
    fit_method = fit_pls(ncomp = 5),
    gh = gh,
    control = ctrl,
    verbose = FALSE,
    seed = 42
  )
}

test_that("plot.mbl rejects non-mbl objects", {
  expect_error(
    plot(structure(list(), class = "not_mbl")),
    "'x' is a list."
  )
})

test_that("plot.mbl validation plot works and restores par", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  obj <- .build_mbl_for_plot(gh = FALSE)
  
  old_par <- par(c("mfrow", "mar"))
  
  .with_test_device({
    expect_invisible(
      plot(obj, what = "validation", metric = "rmse", main = "MBL validation")
    )
  })
  
  expect_equal(par("mfrow"), old_par$mfrow)
  expect_equal(par("mar"), old_par$mar)
})

test_that("plot.mbl GH plot works when gh is available", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  obj <- .build_mbl_for_plot(gh = TRUE)
  
  .with_test_device({
    expect_invisible(
      plot(obj, what = "gh", ncomp = c(1, 2), main = "MBL GH")
    )
  })
})

test_that("plot.mbl combined validation and GH plots work", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  obj <- .build_mbl_for_plot(gh = TRUE)
  
  old_par <- par(c("mfrow", "mar"))
  
  .with_test_device({
    expect_invisible(
      plot(
        obj,
        what = c("validation", "gh"),
        metric = "rmse",
        ncomp = c(1, 2),
        main = "MBL combined"
      )
    )
  })
  
  expect_equal(par("mfrow"), old_par$mfrow)
  expect_equal(par("mar"), old_par$mar)
})

test_that("plot.mbl reports missing GH distances", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  obj <- .build_mbl_for_plot(gh = FALSE)
  
  .with_test_device({
    expect_message(
      plot(obj),
      "GH distance not available in this object."
    )
  })
})

test_that(".plot_gh_1d works on simple scores", {
  xr_scores <- matrix(c(-2, -1, 0, 1), ncol = 1)
  xu_scores <- matrix(c(-1.5, 0.5), ncol = 1)
  
  plot_args <- list(
    pch = 16,
    col.axis = grey(0.3),
    main = "1D GH"
  )
  
  .with_test_device({
    expect_true(
      is.list(
        invisible(
          resemble:::.plot_gh_1d(
            xr_scores = xr_scores,
            xu_scores = xu_scores,
            xr_col = rgb(0, 0, 0.4, 0.5),
            xu_col = rgb(1, 0, 0, 0.5),
            plot_args = plot_args
          )
        )
      )
    )
  })
})

test_that(".plot_gh_2d works on simple scores", {
  xr_scores <- cbind(c(-2, -1, 1, 2), c(-1, 1, -1, 1))
  xu_scores <- cbind(c(-0.5, 0.5), c(0.25, -0.25))
  
  plot_args <- list(
    pch = 16,
    col.axis = grey(0.3),
    main = "2D GH"
  )
  
  .with_test_device({
    expect_null(
      invisible(
        resemble:::.plot_gh_2d(
          xr_scores = xr_scores,
          xu_scores = xu_scores,
          ncomp = c(1, 2),
          xr_col = rgb(0, 0, 0.4, 0.5),
          xu_col = rgb(1, 0, 0, 0.5),
          plot_args = plot_args
        )
      )
    )
  })
})


# =============================================================================
# Deprecated argument tests
# =============================================================================

test_that("mbl rejects deprecated k, k_diss, k_range arguments", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      k = 30,
      verbose = FALSE
    ),
    "Arguments 'k', 'k_diss', 'k_range' have been removed"
  )
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      k_diss = 0.5,
      verbose = FALSE
    ),
    "Arguments 'k', 'k_diss', 'k_range' have been removed"
  )
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      k_range = c(10, 50),
      verbose = FALSE
    ),
    "Arguments 'k', 'k_diss', 'k_range' have been removed"
  )
})

test_that("mbl rejects deprecated method argument", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      method = "pls",
      verbose = FALSE
    ),
    "Argument 'method' has been renamed to 'fit_method'"
  )
})

test_that("mbl rejects deprecated pc_selection argument", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      pc_selection = list(method = "opc"),
      verbose = FALSE
    ),
    "Argument 'pc_selection' has been removed"
  )
})

test_that("mbl rejects deprecated center and scale arguments", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      center = TRUE,
      verbose = FALSE
    ),
    "Arguments 'center' and 'scale' have been removed"
  )
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      scale = FALSE,
      verbose = FALSE
    ),
    "Arguments 'center' and 'scale' have been removed"
  )
})

test_that("mbl rejects deprecated documentation argument", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      documentation = "test",
      verbose = FALSE
    ),
    "Argument 'documentation' has been removed"
  )
})

# =============================================================================
# Input validation tests
# =============================================================================

test_that("mbl requires neighbors argument", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "'neighbors' is required"
  )
})

test_that("mbl validates neighbors class", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = c(20, 30),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "'neighbors' must be created by neighbors_k\\(\\) or neighbors_diss\\(\\)"
  )
})

test_that("mbl validates control class", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      control = list(validation_type = "NNv"),
      verbose = FALSE
    ),
    "'control' must be created by mbl_control\\(\\)"
  )
})

test_that("mbl validates verbose argument", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = "yes"
    ),
    "'verbose' must be TRUE or FALSE"
  )
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = c(TRUE, FALSE)
    ),
    "'verbose' must be TRUE or FALSE"
  )
})

test_that("mbl validates ncol(Xr) == ncol(Xu)", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu[, 1:10],
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "ncol\\(Xr\\) must equal ncol\\(Xu\\)"
  )
})

test_that("mbl validates Xr has at least 4 columns", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  Xr_small <- d$Xr[, 1:3]
  Xu_small <- d$Xu[, 1:3]
  colnames(Xr_small) <- colnames(Xu_small) <- paste0("V", 1:3)
  
  expect_error(
    mbl(
      Xr = Xr_small, Yr = d$Yr, Xu = Xu_small,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 2),
      verbose = FALSE
    ),
    "Xr must have at least 4 columns"
  )
})

test_that("mbl validates Yr length matches Xr rows", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr[1:10], Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "length\\(Yr\\) must equal nrow\\(Xr\\)"
  )
})

test_that("mbl rejects NA values in Yr", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  Yr_na <- d$Yr
  Yr_na[1] <- NA
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = Yr_na, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "NA values in Yr are not supported"
  )
})

test_that("mbl validates Yu length matches Xu rows", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu[1:5],
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "length\\(Yu\\) must equal nrow\\(Xu\\)"
  )
})

test_that("mbl validates diss_method type", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      diss_method = "pca",
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "'diss_method' must be a diss_\\*\\(\\) constructor or a numeric matrix"
  )
})

test_that("mbl validates fit_method type", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = "pls",
      verbose = FALSE
    ),
    "'fit_method' must be created by a fit_\\*\\(\\) constructor"
  )
})

test_that("mbl validates k does not exceed nrow(Xr)", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(100),  # exceeds nrow(Xr) = 40
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "k.*cannot exceed nrow\\(Xr\\)"
  )
})

test_that("mbl validates k_max in neighbors_diss does not exceed nrow(Xr)", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_diss(threshold = 0.5, k_min = 10, k_max = 100),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "k_max.*cannot exceed nrow\\(Xr\\)"
  )
})

test_that("mbl messages when k_max is infinite in neighbors_diss", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_message(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_diss(threshold = 0.5, k_min = 10, k_max = Inf),
      fit_method = fit_pls(ncomp = 5),
      control = mbl_control(validation_type = "NNv", allow_parallel = FALSE),
      verbose = FALSE
    ),
    "setting 'k_max'"
  )
})

test_that("mbl warns when diss_method ncomp exceeds data dimensions", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_warning(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      diss_method = diss_pca(ncomp = 1000),  # exceeds dimensions
      fit_method = fit_pls(ncomp = 5),
      control = mbl_control(validation_type = "NNv", allow_parallel = FALSE),
      verbose = FALSE
    ),
    "Requested ncomp.*exceeds max allowed"
  )
})

test_that("mbl stops when fit_pls ncomp exceeds data dimensions", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 26),  # exceeds dimensions
      control = mbl_control(validation_type = "NNv"),
      verbose = FALSE
    ),
    "More PLS components than observations"
  )
})

test_that("mbl stops when fit_wapls max_ncomp exceeds data dimensions", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 26),  # exceeds dimensions
      control = mbl_control(validation_type = "NNv", allow_parallel = FALSE),
      verbose = FALSE
    ),
    "More PLS components than observations"
  )
})

# =============================================================================
# Precomputed dissimilarity matrix validation
# =============================================================================

test_that("mbl validates precomputed matrix for diss_usage = 'predictors'", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  n_xr <- nrow(d$Xr)
  n_xu <- nrow(d$Xu)
  
  # Wrong dimensions
  bad_diss <- matrix(0, n_xr, n_xu)
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      diss_method = bad_diss,
      diss_usage = "predictors",
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "When diss_usage = 'predictors'.*square"
  )
})

test_that("mbl validates precomputed matrix for diss_usage = 'weights'", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  n_xr <- nrow(d$Xr)
  n_xu <- nrow(d$Xu)
  n_total <- n_xr + n_xu
  
  # Wrong dimensions (square instead of n_xr x n_xu)
  bad_diss <- matrix(0, n_total, n_total)
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      diss_method = bad_diss,
      diss_usage = "weights",
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "When diss_usage = 'weights' or 'none'.*must have nrow\\(Xr\\).*rows"
  )
})

# =============================================================================
# Local CV validation tests
# =============================================================================

test_that("mbl errors when local CV leaves fewer than 3 observations", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  # Very small k with high p leaves too few for hold-out
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(5),  # very small
      fit_method = fit_pls(ncomp = 2),
      control = mbl_control(
        validation_type = "local_cv",
        p = 0.95,  # leaves only ~5% for hold-out
        allow_parallel = FALSE
      ),
      verbose = FALSE
    ),
    "Local cross-validation requires at least 3 observations"
  )
})

test_that("mbl errors when PLS components exceed neighborhood size", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(8),  # small neighborhood
      fit_method = fit_pls(ncomp = 10),  # more components than neighbors
      control = mbl_control(
        validation_type = "local_cv",
        p = 0.5,
        allow_parallel = FALSE
      ),
      verbose = FALSE
    ),
    "More PLS components than observations"
  )
})

# =============================================================================
# return_dissimilarity tests
# =============================================================================

test_that("mbl with return_dissimilarity = TRUE returns dissimilarity", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    return_dissimilarity = TRUE,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
    neighbors = neighbors_k(25),
    diss_method = diss_pca(ncomp = 10),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_true("diss_xr_xu" %in% names(result$dissimilarities))
  expect_true(is.matrix(result$dissimilarities$diss_xr_xu$dissimilarity))
})

test_that("mbl with return_dissimilarity and projection returns projection", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "NNv",
    return_dissimilarity = TRUE,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
    neighbors = neighbors_k(25),
    diss_method = diss_pca(ncomp = 10, return_projection = TRUE),
    diss_usage = "predictors",
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  expect_true("projection" %in% names(result$dissimilarities))
})

# =============================================================================
# get_predictions tests
# =============================================================================

test_that("get_predictions works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
    neighbors = neighbors_k(c(25, 30)),
    fit_method = fit_pls(ncomp = 5),
    control = mbl_control(validation_type = "NNv", allow_parallel = FALSE),
    verbose = FALSE
  )
  
  preds <- get_predictions(result)
  
  expect_true(is.data.frame(preds))
  expect_equal(nrow(preds), nrow(d$Xu))
  expect_equal(ncol(preds), 2)  # Two k values
})

test_that("get_predictions rejects non-mbl objects", {
  expect_error(
    get_predictions(list(a = 1)),
    "the object is not of class 'mbl'"
  )
})

# =============================================================================
# plot.mbl local_cv validation tests
# =============================================================================

test_that("plot.mbl with local_cv validation works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = "local_cv",
    number = 4,
    p = 0.7,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(25),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  .with_test_device({
    expect_invisible(
      plot(result, what = "validation", metric = "rmse")
    )
  })
})

test_that("plot.mbl with both validation types works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  ctrl <- mbl_control(
    validation_type = c("NNv", "local_cv"),
    number = 4,
    p = 0.7,
    allow_parallel = FALSE
  )
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr, Xu = d$Xu, Yu = d$Yu,
    neighbors = neighbors_k(25),
    fit_method = fit_pls(ncomp = 5),
    control = ctrl,
    verbose = FALSE
  )
  
  .with_test_device({
    expect_invisible(
      plot(result, what = "validation", metric = "st_rmse")
    )
  })
})

# =============================================================================
# Column names validation tests
# =============================================================================

test_that("mbl requires Xr to have column names", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  Xr_no_names <- d$Xr
  colnames(Xr_no_names) <- NULL
  
  expect_error(
    mbl(
      Xr = Xr_no_names, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "Xr must have column names"
  )
})

test_that("mbl requires Xu to have column names", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  Xu_no_names <- d$Xu
  colnames(Xu_no_names) <- NULL
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = Xu_no_names,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "Xu must have column names"
  )
})

test_that("mbl requires identical column names in Xr and Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  Xu_diff_names <- d$Xu
  colnames(Xu_diff_names) <- paste0("different_", seq_len(ncol(Xu_diff_names)))
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = Xu_diff_names,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      verbose = FALSE
    ),
    "Column names of Xr and Xu must be identical"
  )
})

# =============================================================================
# group validation tests
# =============================================================================

test_that("mbl validates group length", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu,
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      group = rep(1, 10),  # wrong length
      verbose = FALSE
    ),
    "length\\(group\\) must equal nrow\\(Xr\\)"
  )
})

# =============================================================================
# NNv minimum observations test
# =============================================================================

test_that("mbl with NNv requires at least 3 observations in Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_mbl_data()
  
  expect_error(
    mbl(
      Xr = d$Xr, Yr = d$Yr, Xu = d$Xu[1:2, ],
      neighbors = neighbors_k(25),
      fit_method = fit_pls(ncomp = 5),
      control = mbl_control(validation_type = "NNv", allow_parallel = FALSE),
      verbose = FALSE
    ),
    "Nearest neighbor validation requires at least 3 observations"
  )
})