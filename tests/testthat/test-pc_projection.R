context("test-pc_projection")

test_that("pc_projection works", {
  nirdata <- data("NIRsoil", package = "prospectr")

  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]

  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]

  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]

  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]

  
  cumvar_value <- 0.999
  one_input_matrix <- pc_projection(Xr, 
                                    pc_selection = list(method = "cumvar", value = cumvar_value),
                                    center = TRUE, scale = FALSE,
                                    method = "pca")
  
  expect_true(ncol(one_input_matrix$scores) == one_input_matrix$n_components)
  expect_true(all(one_input_matrix$variance[2,] < cumvar_value))
  
  two_input_matrices <- pc_projection(Xr, Xu,
                                      pc_selection = list(method = "cumvar", value = cumvar_value),
                                      center = TRUE, scale = FALSE,
                                      method = "pca")
  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  expect_true(all(two_input_matrices$variance[2,] < cumvar_value))

  opc_method <- pc_projection(Xr, Xu, Yr = Yr,
                                      pc_selection = list(method = "opc", value = 30),
                                      center = TRUE, scale = FALSE,
                                      method = "pca")
  
  expect_true(opc_method$n_components == which.min(opc_method$opc_evaluation[,2]))
  
  
  expect_true(ncol(two_input_matrices$scores) == two_input_matrices$n_components)
  expect_true(all(two_input_matrices$variance[2,] < cumvar_value))
  
  
})
