context("test-pc_projection")
library(foreach)
library(RhpcBLASctl)
registerDoSEQ()
# tests/testthat/test-ortho_projection.R

# =============================================================================
# Setup helper
# =============================================================================

.setup_ortho_data <- function(n_xr = 40, n_xu = 20) {
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
# Basic ortho_projection tests
# =============================================================================

test_that("ortho_projection with ncomp_by_cumvar works (single matrix)", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  cumvar_value <- 0.999
  
  result <- ortho_projection(
    Xr = d$Xr,
    ncomp = ncomp_by_cumvar(cumvar_value),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
  
  # All but last component should have cumvar < threshold
  test_ncomp <- result$ncomp - 1
  expect_true(all(
    result$variance$x_var["cumulative_explained_var", seq_len(test_ncomp)] < cumvar_value
  ))
})


test_that("ortho_projection with ncomp_by_cumvar works (two matrices)", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  cumvar_value <- 0.999
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    ncomp = ncomp_by_cumvar(cumvar_value),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
  expect_equal(nrow(result$scores), nrow(d$Xr) + nrow(d$Xu))
  
  test_ncomp <- result$ncomp - 1
  expect_true(all(
    result$variance$x_var["cumulative_explained_var", seq_len(test_ncomp)] < cumvar_value
  ))
})


test_that("ortho_projection with ncomp_by_var works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  min_var <- 0.001  # equivalent to 1 - 0.999
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    ncomp = ncomp_by_var(min_var),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(ncol(result$scores), result$ncomp)
  
  # All retained components should explain >= min_var
  expect_true(all(
    result$variance$x_var["explained_var", seq_len(result$ncomp)] >= min_var
  ))
})


test_that("ortho_projection with ncomp_fixed works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    ncomp = ncomp_fixed(10),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(result$ncomp, 10)
  expect_equal(ncol(result$scores), 10)
})


test_that("ortho_projection with integer ncomp works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    ncomp = 10,
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(result$ncomp, 10)
})


# =============================================================================
# predict method tests
# =============================================================================

test_that("predict.ortho_projection works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  tol <- 1e-5
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    ncomp = ncomp_by_cumvar(0.999),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  # Predicting Xr should match stored scores for Xr rows
  preds <- predict(result, d$Xr)
  stored_xr <- predict(result)[seq_len(nrow(d$Xr)), ]
  
  expect_lt(sum(abs(preds - stored_xr)), tol)
})


# =============================================================================
# ncomp_by_opc tests
# =============================================================================

test_that("ortho_projection with ncomp_by_opc works (PCA)", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    center = TRUE, 
    scale = TRUE,
    method = "pca"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_true("opc_evaluation" %in% names(result))
  
  # ncomp should match the minimum RMSD position
  expect_equal(result$ncomp, as.vector(which.min(result$opc_evaluation[, "rmsd_Yr"])))
  expect_equal(result$ncomp, 7)
})


test_that("ortho_projection pca vs pca_nipals give equivalent results", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  tol <- 1e-5
  
  opc_pca <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    center = TRUE, 
    scale = TRUE,
    method = "pca"
  )
  
  opc_nipals <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 30),
    center = TRUE, 
    scale = TRUE,
    method = "pca_nipals"
  )
  
  # Same number of components selected
  expect_equal(opc_pca$ncomp, opc_nipals$ncomp)
  
  # Scores should be equivalent (up to sign)
  cor_equiv <- sapply(
    seq_len(opc_pca$ncomp),
    function(i) abs(cor(opc_pca$scores[, i], opc_nipals$scores[, i]))
  )
  
  expect_lt(sum(1 - cor_equiv), tol)
})


test_that("ncomp_by_var and ncomp_by_cumvar select consistent components", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  cumvar_value <- 0.999
  
  # Get full variance information via OPC (more components)
  opc_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    center = TRUE, 
    scale = TRUE,
    method = "pca"
  )
  
  cumvar_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    ncomp = ncomp_by_cumvar(cumvar_value),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  var_result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    ncomp = ncomp_by_var(1 - cumvar_value),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  # Check consistency between selection methods
  # cumvar: counts components where cumvar < threshold, then adds 1
  expect_equal(
    sum(opc_result$variance$x_var["cumulative_explained_var", ] < cumvar_value),
    cumvar_result$ncomp - 1
  )
  
  # var: counts components where var >= threshold
  expect_equal(
    sum(opc_result$variance$x_var["explained_var", ] >= (1 - cumvar_value)),
    var_result$ncomp
  )
})


# =============================================================================
# Multi-response Yr tests
# =============================================================================

test_that("ortho_projection with multi-column Yr works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  Yr_multi <- cbind(name_test_yr = d$Yr, Yr_2 = d$Yr_2)
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = Yr_multi,
    ncomp = ncomp_by_opc(max_ncomp = 30),
    center = TRUE, 
    scale = FALSE,
    method = "pca_nipals"
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_true("opc_evaluation" %in% names(result))
  expect_true("rmsd_name_test_yr" %in% colnames(result$opc_evaluation))
  expect_true("rmsd_Yr_2" %in% colnames(result$opc_evaluation))
})


# =============================================================================
# PLS method tests
# =============================================================================

test_that("ortho_projection with method = 'pls' works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = ncomp_by_opc(max_ncomp = 15),
    method = "pls",
    scale = TRUE
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_true("weights" %in% names(result))
  expect_true("projection_mat" %in% names(result))
  expect_true("Y_loadings" %in% names(result))
})


test_that("ortho_projection with method = 'mpls' works", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  result <- ortho_projection(
    Xr = d$Xr, 
    Xu = d$Xu,
    Yr = d$Yr,
    ncomp = 10,
    method = "mpls",
    scale = TRUE
  )
  
  expect_s3_class(result, "ortho_projection")
  expect_equal(result$ncomp, 10)
  expect_equal(result$method, "mpls")
})


test_that("ortho_projection PLS requires Yr", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  expect_error(
    ortho_projection(
      Xr = d$Xr, 
      Xu = d$Xu,
      ncomp = 10,
      method = "pls"
    ),
    "'Yr' is required"
  )
})


test_that("ortho_projection ncomp_by_opc requires Yr", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  expect_error(
    ortho_projection(
      Xr = d$Xr, 
      Xu = d$Xu,
      ncomp = ncomp_by_opc(),
      method = "pca"
    ),
    "'Yr' is required"
  )
})


# =============================================================================
# Validation tests
# =============================================================================

test_that("ortho_projection validates center and scale", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  expect_error(
    ortho_projection(d$Xr, center = "yes"),
    "TRUE or FALSE"
  )
  
  expect_error(
    ortho_projection(d$Xr, scale = "yes"),
    "TRUE or FALSE"
  )
})


test_that("ortho_projection validates method", {
  skip_if_not_installed("prospectr")
  
  d <- .setup_ortho_data()
  
  expect_error(
    ortho_projection(d$Xr, method = "invalid"),
    fixed = '"pca", "pca_nipals", "pls", "mpls", "simpls"'
  )
})


# =============================================================================
# Large dataset tests (skipped on CRAN)
# =============================================================================

test_that("ortho_projection works with larger datasets", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  data("NIRsoil", package = "prospectr")
  tol <- 1e-5
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  cumvar_value <- 0.999
  
  # Single matrix
  one_input <- ortho_projection(
    Xr,
    ncomp = ncomp_by_cumvar(cumvar_value),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_equal(ncol(one_input$scores), one_input$ncomp)
  test_ncomp <- one_input$ncomp - 1
  expect_true(all(
    one_input$variance$x_var["cumulative_explained_var", seq_len(test_ncomp)] < cumvar_value
  ))
  
  # Two matrices
  two_input <- ortho_projection(
    Xr, Xu,
    ncomp = ncomp_by_cumvar(cumvar_value),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_equal(ncol(two_input$scores), two_input$ncomp)
  
  # Predictions match
  preds <- sum(abs(
    predict(two_input)[seq_len(nrow(Xr)), ] - predict(two_input, Xr)
  ))
  expect_lt(preds, tol)
  
  # OPC selection
  opc_pca <- ortho_projection(
    Xr, Xu,
    Yr = Yr,
    ncomp = ncomp_by_opc(max_ncomp = 30),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  opc_nipals <- ortho_projection(
    Xr, Xu,
    Yr = Yr,
    ncomp = ncomp_by_opc(max_ncomp = 30),
    center = TRUE, 
    scale = FALSE,
    method = "pca_nipals"
  )
  
  expect_equal(opc_pca$ncomp, as.vector(which.min(opc_pca$opc_evaluation[, "rmsd_Yr"])))
  expect_equal(opc_pca$ncomp, 20)
  
  # PCA and NIPALS should select same number
  expect_equal(opc_pca$ncomp, opc_nipals$ncomp)
  
  # Scores should be equivalent
  cor_equiv <- sapply(
    seq_len(opc_pca$ncomp),
    function(i) abs(cor(opc_pca$scores[, i], opc_nipals$scores[, i]))
  )
  expect_lt(sum(1 - cor_equiv), tol)
  
  # Variance selection consistency
  expect_equal(
    sum(opc_pca$variance$x_var["cumulative_explained_var", ] < cumvar_value),
    two_input$ncomp - 1
  )
  
  var_result <- ortho_projection(
    Xr, Xu,
    ncomp = ncomp_by_var(1 - cumvar_value),
    center = TRUE, 
    scale = FALSE,
    method = "pca"
  )
  
  expect_equal(
    sum(opc_pca$variance$x_var["explained_var", ] >= (1 - cumvar_value)),
    var_result$ncomp
  )
})
