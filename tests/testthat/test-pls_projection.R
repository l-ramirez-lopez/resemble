context("test-pls_projection")

test_that("pls_projection works", {
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

  cumvar_value <- 0.991
  one_input_matrix <- pls_projection(Xr,
    Yr = Yr,
    pc_selection = list(method = "cumvar", value = cumvar_value),
    center = TRUE, scale = FALSE
  )

  test_ncomp <- one_input_matrix$n_components - 1
  expect_true(ncol(one_input_matrix$scores) == one_input_matrix$n_components)
  expect_true(all(one_input_matrix$variance$x_var[3, 1:test_ncomp] < cumvar_value))

  two_input_matrices <- pls_projection(Xr, Xu, Yr,
    pc_selection = list(method = "cumvar", value = cumvar_value),
    center = TRUE, scale = FALSE
  )

  two_input_matrices_var <- pls_projection(Xr, Xu, Yr,
    pc_selection = list(method = "var", value = 1 - cumvar_value),
    center = TRUE, scale = FALSE
  )

  two_test_ncomp <- two_input_matrices$n_components - 1
  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  expect_true(all(two_input_matrices$variance$x_var[3, 1:two_test_ncomp] < cumvar_value))

  preds <- sum(abs(predict(two_input_matrices)[1:nrow(Xr), ] - predict(two_input_matrices, Xr)))
  expect_true(preds < tol)

  opc_method <- pls_projection(Xr, Xu,
    Yr = Yr,
    pc_selection = list(method = "opc", value = 20),
    center = TRUE, scale = FALSE
  )

  expect_true(opc_method$n_components == which.min(opc_method$opc_evaluation[, 2]))
  expect_true(opc_method$n_components == 12)

  # check that the number of components for method = "cumvar" is properly
  # obtained, this can be done with the results of opc_method as it selects more
  # components than in the "cumvar" test
  expect_true(sum(opc_method$variance$x_var[3, ] < cumvar_value) == two_input_matrices$n_components - 1)
  # do the same for method = "var"
  expect_true(sum(opc_method$variance$x_var[2, ] > (1 - cumvar_value)) == two_input_matrices_var$n_components)


  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  two_test_ncnomp <- two_input_matrices$n_components - 1
  expect_true(all(two_input_matrices$variance$x_var[3, 1:two_test_ncnomp] < cumvar_value))

  opc_method <- pls_projection(Xr, Xu,
    Yr = Yr,
    pc_selection = list(method = "opc", value = 20),
    center = TRUE, scale = FALSE
  )

  ## pls2
  pls2_opc_method <- pls_projection(Xr, Xu,
    Yr = cbind(Yr, Yr_2),
    pc_selection = list(method = "opc", value = 20),
    center = TRUE, scale = FALSE
  )

  expect_true(pls2_opc_method$n_components == which.min(pls2_opc_method$opc_evaluation[, 4]))
  distm <- as.matrix(dist(scale(pls2_opc_method$scores[1:nrow(Xr), ], TRUE, TRUE)))
  nn <- apply(distm, MARGIN = 2, FUN = function(x) order(x)[2])

  result_rmsd <- colMeans((cbind(Yr, Yr_2) - cbind(Yr, Yr_2)[nn, ])^2)^0.5
  expect_true(all(pls2_opc_method$opc_evaluation[pls2_opc_method$n_components, 2:3] == result_rmsd))
})
