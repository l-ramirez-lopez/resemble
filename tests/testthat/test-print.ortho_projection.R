context("test-print.ortho_projection")

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
# Tests for print.ortho_projection
# =============================================================================

test_that("print.ortho_projection works for PCA", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, ncomp = 5)
  
  expect_output(print(proj), "Method:")
  expect_output(print(proj), "pca")
  expect_output(print(proj), "Number of components retained:")
  expect_output(print(proj), "5")
  expect_output(print(proj), "Original variance")
  expect_output(print(proj), "Explained variance")
})

test_that("print.ortho_projection works for PCA with Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, Xu = d$Xu, ncomp = 5)
  
  expect_output(print(proj), "Xr; Xu")
})

test_that("print.ortho_projection works for PLS", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, Yr = d$Yr, method = "pls", ncomp = 5)
  
  expect_output(print(proj), "Method:")
  expect_output(print(proj), "pls")
  expect_output(print(proj), "Explained variance in Yr")
})

test_that("print.ortho_projection works for mpls", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, Yr = d$Yr, method = "mpls", ncomp = 5)
  
  expect_output(print(proj), "mpls")
  expect_output(print(proj), "Explained variance in Yr")
})
