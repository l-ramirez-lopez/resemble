context("test-pc_projection")


test_that("pc_projection works", {
  # tolernce for results supposed to be 0s
  tol <- 1e-5
  nirdata <- data("NIRsoil", package = "prospectr")

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

  Xu <- Xu[1:20, ]
  Yu <- Yu[1:20]

  Xr <- Xr[1:40, ]
  Yr <- Yr[1:40]
  Yr_2 <- Yr_2[1:40]

  cumvar_value <- 0.999
  one_input_matrix <- pc_projection(Xr,
    pc_selection = list(method = "cumvar", value = cumvar_value),
    center = TRUE, scale = FALSE,
    method = "pca"
  )

  expect_true(ncol(one_input_matrix$scores) == one_input_matrix$n_components)
  test_ncomp <- one_input_matrix$n_components - 1
  expect_true(all(one_input_matrix$variance$x_var[3, 1:test_ncomp] < cumvar_value))

  two_input_matrices <- pc_projection(Xr, Xu,
    pc_selection = list(method = "cumvar", value = cumvar_value),
    center = TRUE, scale = FALSE,
    method = "pca"
  )

  two_input_matrices_var <- pc_projection(Xr, Xu,
    pc_selection = list(method = "var", value = 1 - cumvar_value),
    center = TRUE, scale = FALSE,
    method = "pca"
  )


  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  two_test_ncomp <- two_input_matrices$n_components - 1
  expect_true(all(two_input_matrices$variance$x_var[3, 1:two_test_ncomp] < cumvar_value))

  preds <- sum(abs(predict(two_input_matrices)[1:nrow(Xr), ] - predict(two_input_matrices, Xr)))
  expect_true(preds < tol)

  opc_method <- pc_projection(Xr, Xu,
    Yr = Yr,
    pc_selection = list(method = "opc", value = 15),
    center = TRUE, scale = TRUE,
    method = "pca"
  )

  opc_method_nipals <- pc_projection(Xr, Xu,
    Yr = Yr,
    pc_selection = list(method = "opc", value = 30),
    center = TRUE, scale = TRUE,
    method = "pca.nipals"
  )

  expect_true(opc_method$n_components == which.min(opc_method$opc_evaluation[, 2]))
  expect_true(opc_method$n_components == 7)

  # check that nipals is equivalent to svd method
  expect_true(opc_method$n_components == opc_method_nipals$n_components)

  cor_equiv <- sapply(1:opc_method$n_components,
    FUN = function(x, y, i) abs(cor(x[, i], y[, i])),
    x = opc_method_nipals$scores,
    y = opc_method$scores
  )

  expect_true(sum(1 - cor_equiv) < tol)

  # check that the number of components for method = "cumvar" is properly
  # obtained, this can be done with the results of opc_method as it selects more
  # components than in the "cumvar" test
  expect_true(sum(opc_method$variance$x_var[3, ] < cumvar_value) == two_input_matrices$n_components - 1)
  # do the same for method = "var"
  expect_true(sum(opc_method$variance$x_var[2, ] > (1 - cumvar_value)) == two_input_matrices_var$n_components)


  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  test_ncomp <- two_input_matrices$n_components - 1
  expect_true(all(two_input_matrices$variance$x_var[3, 1:test_ncomp] < cumvar_value))


  bb <- cbind(name_test_yr = Yr, Yr_2)

  opc_method_nipals <- pc_projection(Xr, Xu,
    Yr = bb,
    pc_selection = list(method = "opc", value = 30),
    center = TRUE, scale = FALSE,
    method = "pca.nipals"
  )

  expect_true("rmsd_name_test_yr" %in% colnames(opc_method_nipals$opc_evaluation))
})


test_that("pc_projection large sets works", {
  skip_on_cran()
  skip_on_travis()
  # tolernce for results supposed to be 0s
  tol <- 1e-5
  nirdata <- data("NIRsoil", package = "prospectr")

  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]

  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Yr_2 <- NIRsoil$Ciso[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]

  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]

  Yu <- Yu[!is.na(Yu)]
  Yr_2 <- Yr[!is.na(Yr)]
  Yr <- Yr[!is.na(Yr)]


  cumvar_value <- 0.999
  one_input_matrix <- pc_projection(Xr,
    pc_selection = list(method = "cumvar", value = cumvar_value),
    center = TRUE, scale = FALSE,
    method = "pca"
  )

  expect_true(ncol(one_input_matrix$scores) == one_input_matrix$n_components)
  test_ncomp <- one_input_matrix$n_components - 1
  expect_true(all(one_input_matrix$variance$x_var[3, 1:test_ncomp] < cumvar_value))

  two_input_matrices <- pc_projection(Xr, Xu,
    pc_selection = list(method = "cumvar", value = cumvar_value),
    center = TRUE, scale = FALSE,
    method = "pca"
  )

  two_input_matrices_var <- pc_projection(Xr, Xu,
    pc_selection = list(method = "var", value = 1 - cumvar_value),
    center = TRUE, scale = FALSE,
    method = "pca"
  )


  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  two_test_ncomp <- two_input_matrices$n_components - 1
  expect_true(all(two_input_matrices$variance$x_var[3, 1:two_test_ncomp] < cumvar_value))

  preds <- sum(abs(predict(two_input_matrices)[1:nrow(Xr), ] - predict(two_input_matrices, Xr)))
  expect_true(preds < tol)

  opc_method <- pc_projection(Xr, Xu,
    Yr = Yr,
    pc_selection = list(method = "opc", value = 30),
    center = TRUE, scale = FALSE,
    method = "pca"
  )

  opc_method_nipals <- pc_projection(Xr, Xu,
    Yr = Yr,
    pc_selection = list(method = "opc", value = 30),
    center = TRUE, scale = FALSE,
    method = "pca.nipals"
  )

  expect_true(opc_method$n_components == which.min(opc_method$opc_evaluation[, 2]))
  expect_true(opc_method$n_components == 20)

  # check that nipals is equivalent to svd method
  expect_true(opc_method$n_components == opc_method_nipals$n_components)

  cor_equiv <- sapply(1:opc_method$n_components,
    FUN = function(x, y, i) abs(cor(x[, i], y[, i])),
    x = opc_method_nipals$scores,
    y = opc_method$scores
  )

  expect_true(sum(1 - cor_equiv) < tol)

  # check that the number of components for method = "cumvar" is properly
  # obtained, this can be done with the results of opc_method as it selects more
  # components than in the "cumvar" test
  expect_true(sum(opc_method$variance$x_var[3, ] < cumvar_value) == two_input_matrices$n_components - 1)
  # do the same for method = "var"
  expect_true(sum(opc_method$variance$x_var[2, ] > (1 - cumvar_value)) == two_input_matrices_var$n_components)


  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  test_ncomp <- two_input_matrices$n_components - 1
  expect_true(all(two_input_matrices$variance$x_var[3, 1:test_ncomp] < cumvar_value))

  bb <- cbind(name_test_yr = Yr, Yr_2)

  opc_method_nipals <- pc_projection(Xr, Xu,
    Yr = bb,
    pc_selection = list(method = "opc", value = 30),
    center = TRUE, scale = FALSE,
    method = "pca.nipals"
  )

  expect_true("rmsd_name_test_yr" %in% colnames(opc_method_nipals$opc_evaluation))
})
