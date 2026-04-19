context("test-liblex")

library(foreach)
library(RhpcBLASctl)
registerDoSEQ()

# =============================================================================
# Setup helper
# =============================================================================

.setup_liblex_data <- function() {
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  
  # Preprocess spectra
  sg_det <- prospectr::savitzkyGolay(
    prospectr::detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1
  test_idx <- NIRsoil$train == 0
  
  list(
    train_x = NIRsoil$spc_pr[train_idx, ],
    train_y = NIRsoil$Ciso[train_idx],
    test_x = NIRsoil$spc_pr[test_idx, ],
    test_y = NIRsoil$Ciso[test_idx]
  )
}

local_liblex_setup <- function(env = parent.frame()) {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  library(prospectr)
  data("NIRsoil", package = "prospectr")
  
  NIRsoil$spc_pr <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  
  env$train_x <- NIRsoil$spc_pr[NIRsoil$train == 1, ]
  env$train_y <- NIRsoil$Ciso[NIRsoil$train == 1]
  env$test_x  <- NIRsoil$spc_pr[NIRsoil$train == 0, ]
  env$test_y  <- NIRsoil$Ciso[NIRsoil$train == 0]
}

check_liblex_predictions <- function(model, test_x, test_y, r2_min = 0.82, rmse_max = 0.6) {
  y_hat <- predict(model, test_x, verbose = FALSE)
  r2 <- cor(y_hat$predictions$pred, test_y, use = "complete.obs")^2
  rmse <- sqrt(mean((y_hat$predictions$pred - test_y)^2, na.rm = TRUE))
  
  expect_true(r2 > r2_min, info = paste("R2 too low:", round(r2, 3)))
  expect_true(rmse < rmse_max, info = paste("RMSE too high:", round(rmse, 3)))
}


# =============================================================================
# Input validation tests - Xr
# =============================================================================

test_that("liblex requires Xr to have column names", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  train_x <- d$train_x
  colnames(train_x) <- NULL
  
  expect_error(
    liblex(
      Xr = train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30)
    ),
    "Xr.*must have column names"
  )
})


test_that("liblex requires Xr to be numeric", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  train_x <- d$train_x
  train_x[1, 1] <- "not_numeric"
  
  expect_error(
    liblex(
      Xr = train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30)
    ),
    "Xr.*must be numeric"
  )
})


# =============================================================================
# Input validation tests - Yr
# =============================================================================

test_that("liblex requires Yr to be numeric", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = as.character(d$train_y),
      neighbors = neighbors_k(30)
    ),
    "Yr.*must be.*numeric"
  )
})


test_that("liblex requires Yr length to match Xr rows", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y[1:10],
      neighbors = neighbors_k(30)
    ),
    "'Yr' must have the same number of observations as 'Xr'"
  )
})


test_that("liblex rejects multi-column Yr", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = cbind(d$train_y, d$train_y),
      neighbors = neighbors_k(30)
    ),
    "'Yr' must be a numeric vector or single-column matrix"
  )
})


# =============================================================================
# Input validation tests - neighbors
# =============================================================================

test_that("liblex requires neighbors argument", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y
    ),
    "neighbors.*required"
  )
})

test_that("liblex with neighbors_diss works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  result <- liblex(
    Xr = d$train_x,
    Yr = d$train_y,
    neighbors = neighbors_diss(
      threshold = seq(0.05, 0.3, length.out = 3),
      k_min = 20,
      k_max = 40
    ),
    diss_method = diss_correlation(ws = 27, scale = TRUE),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 8, method = "mpls"),
    control = liblex_control(tune = TRUE),
    verbose = FALSE
  )
  
  expect_s3_class(result, "liblex")
  expect_true("diss_threshold" %in% names(result$results))
})


test_that("liblex with neighbors_diss predictions work", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  lib <- liblex(
    Xr = d$train_x,
    Yr = d$train_y,
    neighbors = neighbors_diss(
      threshold = c(0.1, 0.2, 0.3),
      k_min = 20,
      k_max = 40
    ),
    diss_method = diss_correlation(ws = 27, scale = TRUE),
    fit_method = fit_pls(ncomp = 6),
    control = liblex_control(tune = FALSE),
    verbose = FALSE
  )
  
  preds <- predict(lib, d$test_x, verbose = FALSE)
  
  expect_true(is.data.frame(preds$predictions))
  expect_equal(nrow(preds$predictions), nrow(d$test_x))
})

test_that("liblex requires neighbors_k object", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = 30
    ),
    "neighbors.*must be.*neighbors_k"
  )
})


test_that("liblex requires neighbors values >= 4", {
  
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(3)
    ),
    "All values in 'k' must be at least 4"
  )
})


# =============================================================================
# Input validation tests - diss_method
# =============================================================================

test_that("liblex validates diss_method type", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      diss_method = "euclidean"
    ),
    "diss_method.*must be.*diss_\\*"
  )
})


test_that("liblex validates precomputed diss_method matrix is numeric", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  n <- nrow(d$train_x)
  diss_mat <- matrix("a", n, n)
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      diss_method = diss_mat
    ),
    "diss_method.*matrix must be numeric"
  )
})


test_that("liblex validates precomputed diss_method matrix is square", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  n <- nrow(d$train_x)
  diss_mat <- matrix(0, n, n - 5)
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      diss_method = diss_mat
    ),
    "diss_method.*matrix must be square"
  )
})


test_that("liblex validates precomputed diss_method matrix dimensions", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  n <- nrow(d$train_x)
  diss_mat <- matrix(0, n - 10, n - 10)
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      diss_method = diss_mat
    ),
    "diss_method.*dimensions must match"
  )
})


# =============================================================================
# Input validation tests - fit_method
# =============================================================================

test_that("liblex validates fit_method type", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      fit_method = "wapls"
    ),
    "fit_method.*must be.*fit_\\*"
  )
})


# =============================================================================
# Input validation tests - anchor_indices
# =============================================================================

test_that("liblex validates anchor_indices are numeric without NA", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      anchor_indices = c(1, 2, NA)
    ),
    "anchor_indices.*numeric vector without NA"
  )
})


test_that("liblex validates anchor_indices are within bounds", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      anchor_indices = c(1, 2, nrow(d$train_x) + 10)
    ),
    "anchor_indices.*must be between 1 and nrow"
  )
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      anchor_indices = c(0, 1, 2)
    ),
    "anchor_indices.*must be between 1 and nrow"
  )
})


test_that("liblex validates anchor_indices does not exceed 90% of data", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  n <- nrow(d$train_x)
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      anchor_indices = seq_len(n)  # 100% of data
    ),
    "anchor_indices.*exceeds 90%"
  )
})


# =============================================================================
# Input validation tests - gh
# =============================================================================

test_that("liblex validates gh parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      gh = "yes"
    ),
    "gh.*must be TRUE or FALSE"
  )
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      gh = NA
    ),
    "gh.*must be TRUE or FALSE"
  )
})


# =============================================================================
# Input validation tests - control
# =============================================================================

test_that("liblex validates control parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      control = list(mode = "build")
    ),
    "control.*must be created by liblex_control"
  )
})


# =============================================================================
# Input validation tests - verbose
# =============================================================================

test_that("liblex validates verbose parameter", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      verbose = "yes"
    ),
    "'verbose' must be a single TRUE or FALSE value"
  )
})


# =============================================================================
# Input validation tests - group
# =============================================================================

test_that("liblex validates group length", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      group = factor(rep(1:5, 10)),
      verbose = FALSE
    ),
    "group.*must have length equal to nrow"
  )
})


# =============================================================================
# Basic functionality tests (skip on CRAN)
# =============================================================================

test_that("liblex runs with default parameters", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  # Use subset for speed
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
  expect_true("dissimilarity" %in% names(model))
  expect_true("fit_method" %in% names(model))
  expect_true("coefficients" %in% names(model))
  expect_true("optimal_params" %in% names(model))
})


test_that("liblex runs with diss_pca", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_pca(ncomp = ncomp_by_var(0.99)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})


test_that("liblex runs with diss_pls", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_pls(ncomp = 10),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})


test_that("liblex runs with diss_correlation", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_correlation(ws = 27, scale = TRUE),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})


test_that("liblex runs with diss_euclidean", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_euclidean(center = TRUE, scale = TRUE),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})


test_that("liblex runs with fit_pls", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_pls(ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})


test_that("liblex runs with anchor_indices", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  anchor_idx <- sample(seq_along(idx), 30)
  
  expect_message(
    model <- liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(c(20, 30)),
      anchor_indices = anchor_idx,
      fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
      verbose = FALSE
    ),
    NA
  )
  
  expect_s3_class(model, "liblex")
  expect_equal(model$anchor_indices, anchor_idx)
})


test_that("liblex runs with gh = FALSE", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    gh = FALSE,
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
  expect_null(model$gh)
})


test_that("liblex runs with gh = TRUE", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    gh = TRUE,
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
  expect_true("gh" %in% names(model))
  expect_true("gh_Xr" %in% names(model$gh))
  expect_true("projection" %in% names(model$gh))
})


test_that("liblex allows missing values in Yr", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  # Introduce some NAs
  train_y <- d$train_y[idx]
  train_y[sample(length(train_y), 5)] <- NA
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = train_y,
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})


test_that("liblex with simpls method works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10, method = "simpls"),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})


# =============================================================================
# Output structure tests
# =============================================================================

test_that("liblex output contains expected elements", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expected_elements <- c(
    "dissimilarity", "fit_method", "gh", "results", "best",
    "optimal_params", "residuals", "coefficients", "vips",
    "selectivity_ratios", "scaling", "neighborhood_stats", "anchor_indices"
  )
  
  for (elem in expected_elements) {
    expect_true(elem %in% names(model), info = paste("Missing element:", elem))
  }
  
  # Check coefficients structure
  expect_true("B0" %in% names(model$coefficients))
  expect_true("B" %in% names(model$coefficients))
  
  # Check optimal_params structure
  expect_true("k" %in% names(model$optimal_params))
})


# =============================================================================
# predict.liblex input validation tests
# =============================================================================


test_that("predict.liblex validates probs parameter", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], probs = c(-0.1, 0.5)),
    "probs.*values in \\[0, 1\\]"
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], probs = c(0.5, 1.5)),
    "probs.*values in \\[0, 1\\]"
  )
})


test_that("predict.liblex validates range_prediction_limits", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], range_prediction_limits = "yes"),
    "range_prediction_limits.*must be TRUE or FALSE"
  )
})


test_that("predict.liblex validates residual_cutoff", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], residual_cutoff = -1),
    "residual_cutoff.*positive"
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], residual_cutoff = 0),
    "residual_cutoff.*positive"
  )
})


test_that("predict.liblex validates enforce_indices", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], enforce_indices = c(0, 1)),
    "enforce_indices.*positive integers"
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], enforce_indices = c(1, 1000)),
    "enforce_indices.*exceeding number of models"
  )
})


test_that("predict.liblex validates verbose", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], verbose = "yes"),
    "verbose.*must be TRUE or FALSE"
  )
})


test_that("predict.liblex validates newdata has required variables", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  # Remove some columns from test data
  test_subset <- d$test_x[1:5, 1:10]
  
  expect_error(
    predict(model, newdata = test_subset, verbose = FALSE),
    "Missing predictor variables"
  )
})


test_that("predict.liblex validates adaptive_bandwidth", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], adaptive_bandwidth = "yes"),
    "adaptive_bandwidth.*must be TRUE or FALSE"
  )
})


test_that("predict.liblex validates reliability_weighting", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], reliability_weighting = "yes"),
    "reliability_weighting.*must be TRUE or FALSE"
  )
})


# =============================================================================
# predict.liblex functionality tests
# =============================================================================

test_that("predict.liblex works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  preds <- predict(model, newdata = d$test_x[1:10, ], verbose = FALSE)
  
  expect_type(preds, "list")
  expect_true("predictions" %in% names(preds))
  expect_true("neighbors" %in% names(preds))
  expect_true("expert_predictions" %in% names(preds))
  
  # Check predictions structure
  expect_true("pred" %in% names(preds$predictions))
  expect_true("pred_sd" %in% names(preds$predictions))
  expect_equal(nrow(preds$predictions), 10)
})


test_that("predict.liblex works with different weighting functions", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  weighting_methods <- c(
    "gaussian", "tricube", "triweight", "triangular",
    "quartic", "parabolic", "cauchy", "none"
  )
  
  for (w in weighting_methods) {
    preds <- predict(
      model, 
      newdata = d$test_x[1:5, ], 
      weighting = w,
      verbose = FALSE
    )
    expect_equal(nrow(preds$predictions), 5, info = paste("Weighting:", w))
  }
})


test_that("predict.liblex includes GH distance when available", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    gh = TRUE,
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  preds <- predict(model, newdata = d$test_x[1:10, ], verbose = FALSE)
  
  expect_true("gh" %in% names(preds$predictions))
  expect_equal(length(preds$predictions$gh), 10)
})


test_that("predict.liblex with range_prediction_limits clips predictions", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(min_ncomp = 3, max_ncomp = 10),
    verbose = FALSE
  )
  
  preds <- predict(
    model, 
    newdata = d$test_x[1:10, ], 
    range_prediction_limits = TRUE,
    verbose = FALSE
  )
  
  # Check that below_min and above_max flags exist
  expect_true("below_min" %in% names(preds$predictions))
  expect_true("above_max" %in% names(preds$predictions))
  expect_true("min_yr" %in% names(preds$predictions))
  expect_true("max_yr" %in% names(preds$predictions))
})



































test_that("liblex with mpls, diss_scale=TRUE, fit_scale=FALSE", {
  skip_on_cran()
  local_liblex_setup()
  
  model <- liblex(
    Xr = train_x, Yr = train_y,
    neighbors = neighbors_k(c(30, 40)),
    diss_method = diss_correlation(ws = 31, scale = TRUE),
    fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, scale = FALSE, method = "mpls"),
    verbose = FALSE,
    control = liblex_control(tune = TRUE)
  )
  
  check_liblex_predictions(model, test_x, test_y)
})

test_that("liblex with simpls, diss_scale=TRUE, fit_scale=FALSE", {
  skip_on_cran()
  local_liblex_setup()
  
  model <- liblex(
    Xr = train_x, Yr = train_y,
    neighbors = neighbors_k(c(30, 40)),
    diss_method = diss_correlation(ws = 31, scale = TRUE),
    fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, scale = FALSE, method = "simpls"),
    verbose = FALSE,
    control = liblex_control(tune = FALSE)
  )
  
  check_liblex_predictions(model, test_x, test_y)
})

test_that("liblex with pls, diss_scale=TRUE, fit_scale=FALSE", {
  skip_on_cran()
  local_liblex_setup()
  
  model <- liblex(
    Xr = train_x, Yr = train_y,
    neighbors = neighbors_k(c(30, 40)),
    diss_method = diss_correlation(ws = 31, scale = TRUE),
    fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, scale = FALSE, method = "pls"),
    verbose = FALSE,
    control = liblex_control(tune = TRUE)
  )
  
  check_liblex_predictions(model, test_x, test_y, r2_min = 0.82, rmse_max = 0.65)
})

test_that("liblex with pls, diss_scale=TRUE, fit_scale=TRUE", {
  skip_on_cran()
  local_liblex_setup()
  
  model <- liblex(
    Xr = train_x, Yr = train_y,
    neighbors = neighbors_k(c(30, 40)),
    diss_method = diss_correlation(ws = 31, scale = TRUE),
    fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, scale = TRUE, method = "pls"),
    verbose = FALSE,
    control = liblex_control(tune = TRUE)
  )

  check_liblex_predictions(model, test_x, test_y, r2_min = 0.82, rmse_max = 0.75)
})

test_that("liblex with pls, diss_scale=FALSE, fit_scale=FALSE", {
  skip_on_cran()
  local_liblex_setup()
  
  model <- liblex(
    Xr = train_x, Yr = train_y,
    neighbors = neighbors_k(c(30, 40)),
    diss_method = diss_correlation(ws = 31, scale = FALSE),
    fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, scale = FALSE, method = "pls"),
    verbose = FALSE,
    control = liblex_control(tune = TRUE)
  )
  check_liblex_predictions(model, test_x, test_y, r2_min = 0.69, rmse_max = 1.3)
})

test_that("liblex with pls, diss_scale=FALSE, fit_scale=TRUE", {
  skip_on_cran()
  local_liblex_setup()
  
  model <- liblex(
    Xr = train_x, Yr = train_y,
    neighbors = neighbors_k(c(30, 40)),
    diss_method = diss_correlation(ws = 31, scale = FALSE),
    fit_method = fit_wapls(min_ncomp = 4, max_ncomp = 17, scale = TRUE, method = "pls"),
    verbose = FALSE,
    control = liblex_control(tune = TRUE)
  )
  check_liblex_predictions(model, test_x, test_y)
})


# =============================================================================
# Additional validation tests for liblex
# =============================================================================

# -----------------------------------------------------------------------------
# anchor_indices validation
# -----------------------------------------------------------------------------

test_that("liblex rejects duplicate anchor_indices", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  
  expect_error(
    liblex(
      Xr = d$train_x,
      Yr = d$train_y,
      neighbors = neighbors_k(30),
      anchor_indices = c(1, 2, 3, 3, 5),
      verbose = FALSE
    ),
    "anchor_indices.*contains duplicate"
  )
})

test_that("liblex warns when max k exceeds anchor count", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  # Use small anchor set with large k
  expect_warning(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(c(20, 25)),
      anchor_indices = 1:15,
      fit_method = fit_wapls(3, 10),
      verbose = FALSE
    ),
    "exceeds number of anchors"
  )
})

# -----------------------------------------------------------------------------
# chunk_size validation
# -----------------------------------------------------------------------------

test_that("liblex errors when chunk_size exceeds data size", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(30),
      fit_method = fit_wapls(3, 10),
      control = liblex_control(chunk_size = 1000),
      verbose = FALSE
    ),
    "chunk_size.*cannot exceed"
  )
})

# -----------------------------------------------------------------------------
# verbose validation (additional coverage)
# -----------------------------------------------------------------------------

test_that("liblex rejects NA verbose", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(30),
      verbose = NA
    ),
    "verbose.*must be "
  )
})

test_that("liblex rejects vector verbose", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(30),
      verbose = c(TRUE, FALSE)
    ),
    "'verbose' must be a single TRUE or FALSE value"
  )
})

# -----------------------------------------------------------------------------
# Precomputed dissimilarity matrix validation
# -----------------------------------------------------------------------------

test_that("liblex validates precomputed matrix row count", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  n <- length(idx)
  
  # Wrong number of rows
  diss_mat <- matrix(0, n - 5, n)
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(30),
      diss_method = diss_mat,
      verbose = FALSE
    ),
    "'diss_method' matrix must be square"
  )
})

test_that("liblex validates precomputed matrix diagonal is zero", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  n <- length(idx)
  
  # Non-zero diagonal
  diss_mat <- matrix(runif(n * n), n, n)
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(30),
      diss_method = diss_mat,
      verbose = FALSE
    ),
    "diagonal.*zeros"
  )
})

test_that("liblex validates precomputed matrix with anchor_indices ncol", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  n <- length(idx)
  anchor_idx <- 1:20
  
  # Wrong number of columns for anchor case
  diss_mat <- matrix(0, n, n)  # Should be n x length(anchor_idx)
  diag(diss_mat) <- 0
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(15),
      anchor_indices = anchor_idx,
      diss_method = diss_mat,
      verbose = FALSE
    ),
    "ncol equal to length\\(anchor_indices\\)"
  )
})

test_that("liblex validates precomputed matrix anchor diagonal is zero", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  n <- length(idx)
  anchor_idx <- 1:20
  
  # Correct dimensions but non-zero anchor diagonal
  diss_mat <- matrix(runif(n * length(anchor_idx)), n, length(anchor_idx))

  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = d$train_y[idx],
      neighbors = neighbors_k(15),
      anchor_indices = anchor_idx,
      diss_method = diss_mat,
      verbose = FALSE
    ),
    "diss_method\\[anchor_indices, \\].*zeros on the diagonal"
  )
})

# -----------------------------------------------------------------------------
# Yr validation for anchors
# -----------------------------------------------------------------------------

test_that("liblex requires at least 3 non-missing Yr for anchors", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  
  # Set almost all Yr to NA
  train_y <- rep(NA_real_, length(idx))
  train_y[1:2] <- d$train_y[idx][1:2]  # Only 2 non-NA
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = train_y,
      neighbors = neighbors_k(30),
      verbose = FALSE
    ),
    "Each 'side_info' variable must have at least 4 non-missing values"
  )
})

# -----------------------------------------------------------------------------
# max_k exceeds valid Yr count
# -----------------------------------------------------------------------------

test_that("liblex errors when max_k exceeds non-missing Yr count", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  
  # Set most Yr to NA so valid count < k
  train_y <- d$train_y[idx]
  train_y[6:50] <- NA
  
  expect_error(
    liblex(
      Xr = d$train_x[idx, ],
      Yr = train_y,
      neighbors = neighbors_k(10),
      verbose = FALSE
    ),
    "no complete element pairs"
  )
})

# -----------------------------------------------------------------------------
# validate mode output structure
# -----------------------------------------------------------------------------

test_that("liblex validate mode returns correct structure", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  result <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(3, 10),
    control = liblex_control(mode = "validate"),
    verbose = FALSE
  )
  
  expect_s3_class(result, "liblex")
  expect_null(result$coefficients)
  expect_true("results" %in% names(result))
  expect_true("best" %in% names(result))
  expect_true("residuals" %in% names(result))
  expect_true("neighborhood_stats" %in% names(result))
})

test_that("liblex validate mode with rownames preserves them", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  train_x <- d$train_x[idx, ]
  rownames(train_x) <- paste0("obs_", seq_len(nrow(train_x)))
  
  result <- liblex(
    Xr = train_x,
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(3, 10),
    control = liblex_control(mode = "validate"),
    verbose = FALSE
  )
  
  # Check that neighborhood_stats has proper row names
  expect_true(all(grepl("^obs_", rownames(result$neighborhood_stats[[1]]))))
})

# =============================================================================
# predict.liblex additional validation tests
# =============================================================================

test_that("predict.liblex errors when precomputed diss requires diss_method", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  n <- length(idx)
  
  # Build with precomputed matrix
  diss_mat <- as.matrix(dist(d$train_x[idx, 1:10]))
  diss_mat <- diss_mat / max(diss_mat)
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(20),
    diss_method = diss_mat,
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], verbose = FALSE),
    "precomputed dissimilarity.*diss_method.*required"
  )
})

test_that("predict.liblex warns when overriding stored diss_method", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_correlation(ws = 27),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_warning(
    predict(
      model, 
      newdata = d$test_x[1:5, ], 
      diss_method = diss_correlation(ws = 31),  # Different
      verbose = FALSE
    ),
    "Overriding stored"
  )
})

test_that("predict.liblex validates diss_method is diss_* object", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:50]
  n <- length(idx)
  
  # Build with precomputed matrix
  diss_mat <- as.matrix(dist(d$train_x[idx, 1:10]))
  diss_mat <- diss_mat / max(diss_mat)
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(20),
    diss_method = diss_mat,
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_error(
    predict(model, newdata = d$test_x[1:5, ], diss_method = "euclidean", verbose = FALSE),
    "diss_method.*must be a diss_\\*"
  )
})

# -----------------------------------------------------------------------------
# verbose output in predict.liblex
# -----------------------------------------------------------------------------

test_that("predict.liblex verbose output for correlation", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_correlation(ws = 27),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_output(
    predict(model, newdata = d$test_x[1:5, ], verbose = TRUE),
    "correlation dissimilarity.*window size"
  )
})

test_that("predict.liblex verbose output for correlation full window", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_correlation(),  # No ws
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_output(
    predict(model, newdata = d$test_x[1:5, ], verbose = TRUE),
    "full window"
  )
})

test_that("predict.liblex verbose output for other diss methods", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    diss_method = diss_pca(ncomp = 10),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_output(
    predict(model, newdata = d$test_x[1:5, ], verbose = TRUE),
    "pca dissimilarity"
  )
})


# -----------------------------------------------------------------------------
# residual_cutoff functionality
# -----------------------------------------------------------------------------

test_that("predict.liblex with residual_cutoff penalizes high-residual models", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  # With residual cutoff
  preds_cutoff <- predict(
    model, 
    newdata = d$test_x[1:5, ], 
    residual_cutoff = 0.1,
    verbose = FALSE
  )
  
  # Without residual cutoff
  preds_no_cutoff <- predict(
    model, 
    newdata = d$test_x[1:5, ], 
    verbose = FALSE
  )
  
  # Predictions may differ due to penalization
  expect_equal(nrow(preds_cutoff$predictions), 5)
  expect_equal(nrow(preds_no_cutoff$predictions), 5)
})

# -----------------------------------------------------------------------------
# enforce_indices functionality
# -----------------------------------------------------------------------------

test_that("predict.liblex with enforce_indices prepends models", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30)),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  enforce_idx <- c(1, 5, 10)
  
  preds <- predict(
    model, 
    newdata = d$test_x[1:5, ],
    enforce_indices = enforce_idx,
    verbose = FALSE
  )
  
  expect_equal(nrow(preds$predictions), 5)
  
  # Check that enforced indices are in neighbors
  idc <- lapply(
    preds$neighbors$indices, 
    FUN = function(x, indices) all(x[1:length(indices)] == indices), 
    indices = enforce_idx
  )
  idc <- do.call(c, idc)
  
  expect_true(all(idc))
})

# =============================================================================
# plot.liblex tests
# =============================================================================

test_that("plot.liblex works with default parameters", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30, 40)),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_silent(plot(model))
})

test_that("plot.liblex works with what = 'rmse'", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30, 40)),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_silent(plot(model, what = "rmse"))
})

test_that("plot.liblex works with what = 'r2'", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30, 40)),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_silent(plot(model, what = "r2"))
})

test_that("plot.liblex works with what = 'residuals'", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 30, 40)),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_silent(plot(model, what = "residuals"))
})

test_that("plot.liblex works with neighbors_diss", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_diss(
      threshold = c(0.1, 0.2, 0.3),
      k_min = 20,
      k_max = 40
    ),
    diss_method = diss_correlation(ws = 27),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_silent(plot(model))
  expect_silent(plot(model, what = "rmse"))
})


test_that("plot.liblex validates what parameter", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20)),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_null(plot(model))

  model$scaling$local_x_center <- model$results <- NULL  # Remove results to trigger error
  expect_error(
    plot(model),
    "Nothing to plot"
  )
})

# =============================================================================
# anchor_indices with non-projection dissimilarity methods
# =============================================================================

test_that("liblex with anchor_indices and diss_euclidean prescales correctly", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  anchor_idx <- sample(seq_along(idx), 30)
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 25)),
    anchor_indices = anchor_idx,
    diss_method = diss_euclidean(center = TRUE, scale = TRUE),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
  expect_equal(model$anchor_indices, anchor_idx)
})

test_that("liblex with anchor_indices and diss_cosine prescales correctly", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_liblex_data()
  idx <- which(!is.na(d$train_y))[1:80]
  anchor_idx <- sample(seq_along(idx), 30)
  
  model <- liblex(
    Xr = d$train_x[idx, ],
    Yr = d$train_y[idx],
    neighbors = neighbors_k(c(20, 25)),
    anchor_indices = anchor_idx,
    diss_method = diss_cosine(center = TRUE, scale = FALSE),
    fit_method = fit_wapls(3, 10),
    verbose = FALSE
  )
  
  expect_s3_class(model, "liblex")
})

