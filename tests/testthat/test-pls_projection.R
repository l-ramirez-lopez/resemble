context("test-pls_projection")
library(foreach)
library(RhpcBLASctl)
registerDoSEQ()
# =============================================================================
# Setup helper
# =============================================================================

.setup_pls_data <- function(n_xr = 40, n_xu = 20) {
  data("NIRsoil", package = "prospectr")
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Yr_2 <- NIRsoil$Ciso[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  y_sel <- !is.na(Yr) & !is.na(Yr_2)
  Xr <- Xr[y_sel, ]
  
  Yu <- Yu[!is.na(Yu)]
  Yr_2 <- Yr_2[y_sel]
  Yr <- Yr[y_sel]
  
  Xu <- Xu[seq_len(n_xu), ]
  Yu <- Yu[seq_len(n_xu)]
  Xr <- Xr[seq_len(n_xr), ]
  Yr <- Yr[seq_len(n_xr)]
  Yr_2 <- Yr_2[seq_len(n_xr)]
  
  list(Xr = Xr, Xu = Xu, Yr = Yr, Yr_2 = Yr_2, Yu = Yu)
}


# =============================================================================
# PLS with multi-response Y tests
# =============================================================================

test_that("ortho_projection PLS with multi-response Y and ncomp_by_var works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  Yr <- NIRsoil[as.logical(NIRsoil$train), c("Ciso", "Nt")]
  
  # Remove rows with NA in either response
  complete <- complete.cases(Yr)
  Xr <- Xr[complete, ]
  Yr <- Yr[complete, ]
  
  result <- ortho_projection(
    Xr = Xr,
    Yr = as.matrix(Yr),
    ncomp = ncomp_by_var(0.01),
    method = "pls"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
})


test_that("ortho_projection PLS with multi-response Y and ncomp_by_cumvar works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  Yr <- NIRsoil[as.logical(NIRsoil$train), c("Ciso", "Nt")]
  
  complete <- complete.cases(Yr)
  Xr <- Xr[complete, ]
  Yr <- Yr[complete, ]
  
  result <- ortho_projection(
    Xr = Xr,
    Yr = as.matrix(Yr),
    ncomp = ncomp_by_cumvar(0.99),
    method = "pls"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
})


# =============================================================================
# Basic PLS projection tests
# =============================================================================

test_that("ortho_projection with method = 'pls' and ncomp_by_cumvar works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()
  cumvar_value <- 0.74
  
  result <- ortho_projection(
    Xr = d$Xr,
    Yr = d$Yr,
    ncomp = ncomp_by_cumvar(cumvar_value),
    method = "pls",
    scale = FALSE
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
  
  # All but last component should have cumvar > threshold (for PLS, variance explained grows)
  test_ncomp <- result$ncomp - 1
  if (test_ncomp > 0) {
    expect_true(all(
      result$variance$x_var["cumulative_explained_var_X", seq_len(test_ncomp)] < cumvar_value
    ))
  }
})


test_that("ortho_projection PLS with two matrices works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()
  cumvar_value <- 0.74
  tol <- 1e-5
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_cumvar(cumvar_value),
    method = "pls",
    scale = FALSE
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
  expect_equal(nrow(result$scores), nrow(d$Xr) + nrow(d$Xu))
  
  # Predictions for Xr should match stored scores
  preds <- sum(abs(
    predict(result)[seq_len(nrow(d$Xr)), ] - predict(result, d$Xr)
  ))
  expect_lt(preds, tol)
})


test_that("ortho_projection PLS with ncomp_by_var works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()
  min_var <- 0.26  # 1 - 0.74
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_var(min_var),
    method = "pls",
    scale = FALSE
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
})


# =============================================================================
# OPC selection tests
# =============================================================================

test_that("ortho_projection PLS with ncomp_by_opc works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    method = "pls",
    scale = TRUE
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_true("opc_evaluation" %in% names(result))
  
  # ncomp should match minimum RMSD position
  expect_equal(result$ncomp, as.vector(which.min(result$opc_evaluation[, "rmsd_Yr"])))
  expect_equal(result$ncomp, 8)
})


test_that("ortho_projection PLS variance selection consistency", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()
  cumvar_value <- 0.74
  
  # Get full variance info via OPC
  opc_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    method = "pls",
    scale = TRUE
  )
  
  cumvar_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_cumvar(cumvar_value),
    method = "pls",
    scale = FALSE
  )
  
  var_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_var(1 - cumvar_value),
    method = "pls",
    scale = FALSE
  )
  
  # Check consistency
  expect_equal(
    sum(opc_result$variance$x_var["cumulative_explained_var_X", ] < cumvar_value),
    cumvar_result$ncomp - 1
  )
  
  expect_equal(
    sum(opc_result$variance$x_var["explained_var_X", ] >= (1 - cumvar_value)),
    var_result$ncomp
  )
})


# =============================================================================
# PLS2 (multi-response) with OPC tests
# =============================================================================

test_that("ortho_projection PLS2 with ncomp_by_opc works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()
  tol <- 1e-6
  Yr_multi <- cbind(d$Yr, d$Yr_2)
  colnames(Yr_multi) <- c("Yr", "Yr_2")
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = Yr_multi,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    method = "pls",
    scale = TRUE
  )
  
  expect_lt(
    sum(abs(result$scores - predict(result, rbind(d$Xr, d$Xu)))), 
    tol,
    label = "proper pls projection of Xr should match stored scores"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_true("opc_evaluation" %in% names(result))
  
  # Should select based on mean standardized RMSD
  expect_equal(
    result$ncomp, 
    as.vector(which.min(result$opc_evaluation[, "mean_standardized_rmsd_Yr"]))
  )
  
  # Verify RMSD calculation manually
  scores_xr <- result$scores[seq_len(nrow(d$Xr)), ]
  distm <- as.matrix(dist(scale(scores_xr, center = TRUE, scale = TRUE)))
  nn <- apply(distm, MARGIN = 2, FUN = function(x) order(x)[2])
  
  expected_rmsd <- round(
    colMeans((Yr_multi - Yr_multi[nn, ])^2)^0.5, 
    4
  )
  
  actual_rmsd <- round(
    result$opc_evaluation[result$ncomp, c("rmsd_Yr", "rmsd_Yr_2")], 
    4
  )
  
  expect_equal(as.vector(actual_rmsd), as.vector(expected_rmsd))
})


# =============================================================================
# SIMPLS tests
# =============================================================================

test_that("ortho_projection with method = 'simpls' works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  tol <- 1e-6
  d <- .setup_pls_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    method = "simpls",
    scale = TRUE
  )
  
  expect_lt(
    sum(abs(result$scores - predict(result, rbind(d$Xr, d$Xu)))), 
    tol,
    label = "proper simpls projection of Xr should match stored scores"
  )
  expect_s3_class(result, "ortho_projection")
  expect_equal(result$method, "simpls")
  expect_true("opc_evaluation" %in% names(result))
})


test_that("ortho_projection PLS and SIMPLS produce equivalent results", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()

  pls_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = 10,
    method = "pls",
    scale = TRUE
  )
  
  simpls_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = 10,
    method = "simpls",
    scale = TRUE
  )
  
  # Predictions should match
  pls_pred <- predict(pls_result, d$Xu)
  simpls_pred <- predict(simpls_result, d$Xu)
  
  for (i in 1:ncol(pls_pred)) {
    expect_gt(
      cor(pls_pred[, i], simpls_pred[, i])^2,
      0.999,
      label = paste("SIMPLS vs NIPALS PLS: projected scores for factor", i, "should corrrelate")
    )
  }
})


# =============================================================================
# Modified PLS tests
# =============================================================================

test_that("ortho_projection with method = 'mpls' works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  d <- .setup_pls_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = 10,
    method = "mpls",
    scale = TRUE
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(result$method, "mpls")
  expect_equal(result$ncomp, 10)
})


# =============================================================================
# Large dataset tests (skipped on CRAN)
# =============================================================================

test_that("ortho_projection PLS works with larger datasets", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  tol <- 1e-5
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Yr_2 <- NIRsoil$Ciso[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  y_sel <- !is.na(Yr) & !is.na(Yr_2)
  Xr <- Xr[y_sel, ]
  Yu <- Yu[!is.na(Yu)]
  Yr_2 <- Yr_2[y_sel]
  Yr <- Yr[y_sel]
  
  cumvar_value <- 0.991
  
  # Single matrix
  one_input <- ortho_projection(
    Xr = Xr,
    Yr = Yr,
    ncomp = ncomp_by_cumvar(cumvar_value),
    method = "simpls",
    scale = FALSE
  )

  expect_equal(ncol(one_input$scores), one_input$ncomp)
  test_ncomp <- one_input$ncomp - 1
  expect_true(all(
    one_input$variance$x_var["cumulative_explained_var_X", seq_len(test_ncomp)] < cumvar_value
  ))
  
  # Two matrices
  two_input <- ortho_projection(
    Xr = Xr, 
    Xu = Xu,
    Yr = Yr,
    ncomp = ncomp_by_cumvar(cumvar_value),
    method = "pls",
    scale = FALSE
  )
  
  expect_equal(ncol(two_input$scores), two_input$ncomp)
  
  # Predictions match
  preds <- sum(abs(
    predict(two_input)[seq_len(nrow(Xr)), ] - predict(two_input, Xr)
  ))
  expect_lt(preds, tol)
  
  # OPC selection
  opc_result <- ortho_projection(
    Xr = Xr, 
    Xu = Xu,
    Yr = Yr,
    ncomp = ncomp_by_opc(max_ncomp = 20),
    method = "pls",
    scale = FALSE
  )
  
  expect_equal(opc_result$ncomp, as.vector(which.min(opc_result$opc_evaluation[, "rmsd_Yr"])))
  expect_equal(opc_result$ncomp, 12)
  
  # Variance selection consistency
  var_result <- ortho_projection(
    Xr = Xr, 
    Xu = Xu,
    Yr = Yr,
    ncomp = ncomp_by_var(1 - cumvar_value),
    method = "pls",
    scale = FALSE
  )
  
  expect_equal(
    sum(opc_result$variance$x_var["cumulative_explained_var_X", ] < cumvar_value),
    two_input$ncomp - 1
  )
  
  expect_equal(
    sum(opc_result$variance$x_var["explained_var_X", ] >= (1 - cumvar_value)),
    var_result$ncomp
  )
  
  # PLS2
  Yr_multi <- cbind(Yr, Yr_2)
  colnames(Yr_multi) <- c("Yr", "Yr_2")
  
  pls2_result <- ortho_projection(
    Xr = Xr, 
    Xu = Xu,
    Yr = Yr_multi,
    ncomp = ncomp_by_opc(max_ncomp = 20),
    method = "pls",
    scale = FALSE
  )
  
  expect_equal(
    pls2_result$ncomp, 
    as.vector(which.min(pls2_result$opc_evaluation[, "mean_standardized_rmsd_Yr"]))
  )
  
  # Verify RMSD calculation
  scores_xr <- pls2_result$scores[seq_len(nrow(Xr)), ]
  distm <- as.matrix(dist(scale(scores_xr, center = TRUE, scale = TRUE)))
  nn <- apply(distm, MARGIN = 2, FUN = function(x) order(x)[2])
  
  expected_rmsd <- round(
    colMeans((Yr_multi - Yr_multi[nn, ])^2)^0.5, 
    4
  )
  
  actual_rmsd <- round(
    pls2_result$opc_evaluation[pls2_result$ncomp, c("rmsd_Yr", "rmsd_Yr_2")], 
    4
  )
  
  expect_equal(as.vector(actual_rmsd), as.vector(expected_rmsd))
})
