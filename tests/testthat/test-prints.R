# =============================================================================
# Tests for print helpers and print methods
# =============================================================================

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


# -----------------------------------------------------------------------------
# Helper function tests
# -----------------------------------------------------------------------------

test_that(".divider returns correct width", {
  skip_on_cran()
  
  div <- .divider(55)
  expect_equal(nchar(div), 55)
  expect_true(all(strsplit(div, "")[[1]] == "_"))
  
  div_short <- .divider(10)
  expect_equal(nchar(div_short), 10)
})

test_that(".divider respects system width", {
  skip_on_cran()
  
  old_width <- getOption("width")
  on.exit(options(width = old_width))
  
  options(width = 40)
  div <- .divider(100)
  expect_equal(nchar(div), 40)
})

test_that(".use_color returns logical", {
  skip_on_cran()
  
  result <- .use_color()
  expect_type(result, "logical")
})

test_that(".col_blue returns string", {
  skip_on_cran()
  
  result <- .col_blue("test")
  expect_type(result, "character")
  expect_true(grepl("test", result))
})

test_that(".col_bold_red returns string", {
  skip_on_cran()
  
  result <- .col_bold_red("test")
  expect_type(result, "character")
  expect_true(grepl("test", result))
})

test_that(".truncate_call truncates long calls", {
  skip_on_cran()
  
  short_call <- quote(foo(x = 1))
  result_short <- .truncate_call(short_call)
  expect_type(result_short, "character")
  
  long_call <- quote(foo(a = 1, b = 2, c = 3, d = 4, e = 5, f = 6, g = 7, h = 8))
  result_long <- .truncate_call(long_call, width = 20)
  expect_true(length(result_long) <= 4)
})

# -----------------------------------------------------------------------------
# print.ortho_projection tests
# -----------------------------------------------------------------------------

test_that("print.ortho_projection works for PCA", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, ncomp = 5)
  out <- capture.output(print(proj))
  
  expect_true(any(grepl("Method:", out)))
  expect_true(any(grepl("pca", out)))
  expect_true(any(grepl("Number of components retained:", out)))
  expect_true(any(grepl("Original variance", out)))
  expect_true(any(grepl("Explained variance", out)))
})

test_that("print.ortho_projection works for PCA with Xu", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, Xu = d$Xu, ncomp = 5)
  out <- capture.output(print(proj))
  
  expect_true(any(grepl("Xr; Xu", out)))
})

test_that("print.ortho_projection works for PLS", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, Yr = d$Yr, method = "pls", ncomp = 5)
  out <- capture.output(print(proj))
  
  expect_true(any(grepl("pls", out)))
  expect_true(any(grepl("Explained variance in Yr", out)))
})

test_that("print.ortho_projection works for mpls", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, Yr = d$Yr, method = "mpls", ncomp = 5)
  out <- capture.output(print(proj))
  
  expect_true(any(grepl("mpls", out)))
  expect_true(any(grepl("Explained variance in Yr", out)))
})

test_that("print.ortho_projection returns invisible", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  proj <- ortho_projection(d$Xr, ncomp = 5)
  
  capture.output(expect_invisible(print(proj)))
})

# -----------------------------------------------------------------------------
# print.liblex tests
# -----------------------------------------------------------------------------

test_that("print.liblex works for model library with wapls", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  lib <- liblex(
    Xr = d$Xr, Yr = d$Yr,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(k = c(20, 30)),
    fit_method = fit_wapls(3, 10),
    mode = "build",
    verbose = FALSE
  )
  out <- capture.output(print(lib))
  
  expect_true(any(grepl("liblex model library", out)))
  expect_true(any(grepl("Models:", out)))
  expect_true(any(grepl("Dissimilarity", out)))
  expect_true(any(grepl("Local fit method", out)))
  expect_true(any(grepl("Optimal parameters", out)))
})

test_that("print.liblex works for model library with pls", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  lib <- liblex(
    Xr = d$Xr, Yr = d$Yr,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(k = c(20, 30)),
    fit_method = fit_pls(ncomp = 10),
    mode = "build",
    verbose = FALSE
  )
  out <- capture.output(print(lib))
  
  expect_true(any(grepl("liblex model library", out)))
  expect_true(any(grepl("Models:", out)))
  expect_true(any(grepl("Predictors:", out)))
  expect_true(any(grepl("Dissimilarity", out)))
  expect_true(any(grepl("Local fit method", out)))
  expect_true(any(grepl("Optimal parameters", out)))
})

test_that("print.liblex returns invisible", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  lib <- liblex(
    Xr = d$Xr, Yr = d$Yr,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(k = c(20, 30)),
    fit_method = fit_pls(ncomp = 10),
    mode = "build",
    verbose = FALSE
  )
  
  capture.output(expect_invisible(print(lib)))
})

# -----------------------------------------------------------------------------
# print.mbl tests
# -----------------------------------------------------------------------------

test_that("print.mbl works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr,
    Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(k = 30),
    fit_method = fit_pls(ncomp = 10),
    control = mbl_control(validation_type = "NNv"),
    verbose = FALSE
  )
  out <- capture.output(print(result))
  
  expect_true(any(grepl("mbl predictions", out)))
  expect_true(any(grepl("Predictions:", out)))
  expect_true(any(grepl("Dissimilarity", out)))
  expect_true(any(grepl("Local fit method", out)))
})

test_that("print.mbl returns invisible", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- mbl(
    Xr = d$Xr, Yr = d$Yr,
    Xu = d$Xu,
    diss_method = diss_pca(ncomp = 10),
    neighbors = neighbors_k(k = 30),
    fit_method = fit_pls(ncomp = 10),
    control = mbl_control(validation_type = "NNv"),
    verbose = FALSE
  )
  
  capture.output(expect_invisible(print(result)))
})

# -----------------------------------------------------------------------------
# print.resemble_model tests
# -----------------------------------------------------------------------------

test_that("print.resemble_model works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  mod <- model(
    Xr = d$Xr, Yr = d$Yr,
    fit_method = fit_pls(ncomp = 10),
    verbose = FALSE
  )
  out <- capture.output(print(mod))
  
  expect_true(any(grepl("Global resemble model", out)))
  expect_true(any(grepl("Method:", out)))
  expect_true(any(grepl("Observations:", out)))
  expect_true(any(grepl("Variables:", out)))
  expect_true(any(grepl("Fit method", out)))
  expect_true(any(grepl("Cross-validation", out)))
})

test_that("print.resemble_model returns invisible", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  mod <- model(
    Xr = d$Xr, Yr = d$Yr,
    fit_method = fit_pls(ncomp = 10),
    verbose = FALSE
  )
  
  capture.output(expect_invisible(print(mod)))
})

# -----------------------------------------------------------------------------
# print.gesearch tests
# -----------------------------------------------------------------------------

test_that("print.gesearch works", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- gesearch(
    Xr = d$Xr, Yr = d$Yr,
    Xu = d$Xu, Yu = d$Yu,
    k = 20, b = 20,
    fit_method = fit_pls(ncomp = 10),
    verbose = FALSE
  )
  out <- capture.output(print(result))
  
  expect_true(any(grepl("gesearch results", out)))
  expect_true(any(grepl("Iterations:", out)))
  expect_true(any(grepl("Selected:", out)))
  expect_true(any(grepl("Removed:", out)))
  expect_true(any(grepl("Fit method", out)))
})

test_that("print.gesearch returns invisible", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  d <- .setup_nirsoil_data()
  
  result <- gesearch(
    Xr = d$Xr, Yr = d$Yr,
    Xu = d$Xu, Yu = d$Yu,
    k = 30, b = 10,
    fit_method = fit_pls(ncomp = 10),
    verbose = FALSE
  )
  
  capture.output(expect_invisible(print(result)))
})
