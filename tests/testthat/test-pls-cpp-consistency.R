context("test-pls-cpp-consistency")

test_that("PLS (NIPALS) implementations produce consistent results with scale = TRUE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40
  tol <- 1e-8
  tol2 <- 1e-5
  
  # Fit models
  pls_basics <- resemble:::opls_get_basics(
    X = train_x, 
    Y = as.matrix(train_y), 
    ncomp = ncomp_test,
    scale = TRUE,
    maxiter = 100,
    tol = 1e-06, 
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_opls <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_projection <- resemble:::opls_for_projection(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,            
    maxiter = 1000,
    tol = 1e-06,
    pcSelmethod = "manual",
    pcSelvalue = ncomp_test, 
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_gesearch <- resemble:::opls_gesearch(
    Xr = train_x, 
    Yr = as.matrix(train_y), 
    Xu = test_x, 
    ncomp = ncomp_test,
    scale = TRUE,     
    response = TRUE, 
    reconstruction = TRUE,
    similarity = TRUE,
    fresponse = TRUE,
    algorithm = "pls"
  )
  
  # Predictions
  pls_basics_predict <- resemble:::predict_opls(
    bo = pls_basics$bo, 
    b = pls_basics$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = pls_basics$transf$Xscale
  )
  
  pls_opls_predict <- resemble:::predict_opls(
    bo = pls_opls$bo, 
    b = pls_opls$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = pls_opls$transf$Xscale
  )
  
  # Reconstruction error
  pls_opls_rec <- resemble:::reconstruction_error(
    x = test_x, 
    projection_mat = pls_opls$projection_mat, 
    xloadings = pls_opls$X_loadings,
    scale = TRUE,
    Xcenter = pls_opls$transf$Xcenter,
    Xscale = pls_opls$transf$Xscale,
    scale_back = TRUE
  )
  
  # Tests
  expect_lt(
    sum(abs(pls_projection$weights - pls_basics$weights)), 
    tol,
    label = "PLS scale=TRUE: opls_for_projection vs opls_get_basics weights"
  )
  
  expect_lt(
    sum(abs(pls_basics$weights - pls_opls$weights)), 
    tol,
    label = "PLS scale=TRUE: opls_get_basics vs opls weights"
  )
  
  expect_lt(
    sum(abs(pls_basics_predict - pls_opls_predict)), 
    tol,
    label = "PLS scale=TRUE: predict_opls consistency"
  )
  
  expect_lt(
    sum(abs(pls_projection$scores - pls_opls$scores)), 
    tol,
    label = "PLS scale=TRUE: opls_for_projection vs opls scores"
  )
  
  expect_lt(
    sum(abs(pls_gesearch$pred_response - pls_opls_predict[, ncomp_test])), 
    tol,
    label = "PLS scale=TRUE: opls_gesearch vs opls predictions"
  )
  
  expect_lt(
    abs(drop(pls_gesearch$rmse_reconstruction) - drop(pls_opls_rec)), 
    tol,
    label = "PLS scale=TRUE: opls_gesearch vs reconstruction_error"
  )
})


test_that("PLS (NIPALS) implementations produce consistent results with scale = FALSE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40
  tol <- 1e-8
  
  # Fit models
  pls_basics <- resemble:::opls_get_basics(
    X = train_x, 
    Y = as.matrix(train_y), 
    ncomp = ncomp_test,
    scale = FALSE,
    maxiter = 100,
    tol = 1e-06, 
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_opls <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_projection <- resemble:::opls_for_projection(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,            
    maxiter = 1000,
    tol = 1e-06,
    pcSelmethod = "manual",
    pcSelvalue = ncomp_test, 
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_gesearch <- resemble:::opls_gesearch(
    Xr = train_x, 
    Yr = as.matrix(train_y), 
    Xu = test_x, 
    ncomp = ncomp_test,
    scale = FALSE,     
    response = TRUE, 
    reconstruction = TRUE,
    similarity = TRUE,
    fresponse = TRUE,
    algorithm = "pls"
  )
  
  # Predictions
  pls_basics_predict <- resemble:::predict_opls(
    bo = pls_basics$bo, 
    b = pls_basics$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = pls_basics$transf$Xscale
  )
  
  pls_opls_predict <- resemble:::predict_opls(
    bo = pls_opls$bo, 
    b = pls_opls$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = pls_opls$transf$Xscale
  )
  
  # Reconstruction error
  pls_opls_rec <- resemble:::reconstruction_error(
    x = test_x, 
    projection_mat = pls_opls$projection_mat, 
    xloadings = pls_opls$X_loadings,
    scale = FALSE,
    Xcenter = pls_opls$transf$Xcenter,
    Xscale = pls_opls$transf$Xscale,
    scale_back = TRUE
  )
  
  # Tests
  expect_lt(
    sum(abs(pls_projection$weights - pls_basics$weights)), 
    tol,
    label = "PLS scale=FALSE: opls_for_projection vs opls_get_basics weights"
  )
  
  expect_lt(
    sum(abs(pls_basics$weights - pls_opls$weights)), 
    tol,
    label = "PLS scale=FALSE: opls_get_basics vs opls weights"
  )
  
  expect_lt(
    sum(abs(pls_basics_predict - pls_opls_predict)), 
    tol,
    label = "PLS scale=FALSE: predict_opls consistency"
  )
  
  expect_lt(
    sum(abs(pls_projection$scores - pls_opls$scores)), 
    tol,
    label = "PLS scale=FALSE: opls_for_projection vs opls scores"
  )
  
  expect_lt(
    sum(abs(pls_gesearch$pred_response - pls_opls_predict[, ncomp_test])), 
    tol,
    label = "PLS scale=FALSE: opls_gesearch vs opls predictions"
  )
  
  expect_lt(
    abs(drop(pls_gesearch$rmse_reconstruction) - drop(pls_opls_rec)), 
    tol,
    label = "PLS scale=FALSE: opls_gesearch vs reconstruction_error"
  )
  
  # Xscale should be all ones when scale = FALSE
  expect_true(
    all(pls_basics$transf$Xscale == 1),
    label = "PLS scale=FALSE: Xscale should be all ones (opls_get_basics)"
  )
  
  expect_true(
    all(pls_opls$transf$Xscale == 1),
    label = "PLS scale=FALSE: Xscale should be all ones (opls)"
  )
})


test_that("SIMPLS implementations produce consistent results with scale = TRUE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40
  tol <- 1e-8
  tol2 <- 1e-5
  
  # Fit models
  simpls_basics <- resemble:::opls_get_basics(
    X = train_x, 
    Y = as.matrix(train_y), 
    ncomp = ncomp_test,
    scale = TRUE,
    maxiter = 100,
    tol = 1e-06, 
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_opls <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_projection <- resemble:::opls_for_projection(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,            
    maxiter = 1000,
    tol = 1e-06,
    pcSelmethod = "manual",
    pcSelvalue = ncomp_test, 
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_gesearch <- resemble:::opls_gesearch(
    Xr = train_x, 
    Yr = as.matrix(train_y), 
    Xu = test_x, 
    ncomp = ncomp_test,
    scale = TRUE,     
    response = TRUE, 
    reconstruction = TRUE,
    similarity = TRUE,
    fresponse = TRUE,
    algorithm = "simpls"
  )
  
  # Predictions
  simpls_basics_predict <- resemble:::predict_opls(
    bo = simpls_basics$bo, 
    b = simpls_basics$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = simpls_basics$transf$Xscale
  )
  
  simpls_opls_predict <- resemble:::predict_opls(
    bo = simpls_opls$bo, 
    b = simpls_opls$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = simpls_opls$transf$Xscale
  )
  
  # Reconstruction error
  simpls_opls_rec <- resemble:::reconstruction_error(
    x = test_x, 
    projection_mat = simpls_opls$projection_mat, 
    xloadings = simpls_opls$X_loadings,
    scale = TRUE,
    Xcenter = simpls_opls$transf$Xcenter,
    Xscale = simpls_opls$transf$Xscale,
    scale_back = TRUE
  )
  
  # Tests
  expect_lt(
    sum(abs(simpls_projection$weights - simpls_basics$weights)), 
    tol,
    label = "SIMPLS scale=TRUE: opls_for_projection vs opls_get_basics weights"
  )
  
  expect_lt(
    sum(abs(simpls_basics$weights - simpls_opls$weights)), 
    tol,
    label = "SIMPLS scale=TRUE: opls_get_basics vs opls weights"
  )
  
  expect_lt(
    sum(abs(simpls_basics_predict - simpls_opls_predict)), 
    tol2,
    label = "SIMPLS scale=TRUE: predict_opls consistency"
  )
  
  expect_lt(
    sum(abs(simpls_projection$scores - simpls_opls$scores)), 
    tol,
    label = "SIMPLS scale=TRUE: opls_for_projection vs opls scores"
  )
  
  expect_lt(
    sum(abs(simpls_gesearch$pred_response - simpls_opls_predict[, ncomp_test])), 
    tol2,
    label = "SIMPLS scale=TRUE: opls_gesearch vs opls predictions"
  )
  
  expect_lt(
    abs(drop(simpls_gesearch$rmse_reconstruction) - drop(simpls_opls_rec)), 
    tol,
    label = "SIMPLS scale=TRUE: opls_gesearch vs reconstruction_error"
  )
})


test_that("SIMPLS implementations produce consistent results with scale = FALSE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40
  tol <- 1e-8
  tol2 <- 1e-5
  
  # Fit models
  simpls_basics <- resemble:::opls_get_basics(
    X = train_x, 
    Y = as.matrix(train_y), 
    ncomp = ncomp_test,
    scale = FALSE,
    maxiter = 100,
    tol = 1e-06, 
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_opls <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_projection <- resemble:::opls_for_projection(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,            
    maxiter = 1000,
    tol = 1e-06,
    pcSelmethod = "manual",
    pcSelvalue = ncomp_test, 
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_gesearch <- resemble:::opls_gesearch(
    Xr = train_x, 
    Yr = as.matrix(train_y), 
    Xu = test_x, 
    ncomp = ncomp_test,
    scale = FALSE,     
    response = TRUE, 
    reconstruction = TRUE,
    similarity = TRUE,
    fresponse = TRUE,
    algorithm = "simpls"
  )
  
  # Predictions
  simpls_basics_predict <- resemble:::predict_opls(
    bo = simpls_basics$bo, 
    b = simpls_basics$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = simpls_basics$transf$Xscale
  )
  
  simpls_opls_predict <- resemble:::predict_opls(
    bo = simpls_opls$bo, 
    b = simpls_opls$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = simpls_opls$transf$Xscale
  )
  
  # Reconstruction error
  simpls_opls_rec <- resemble:::reconstruction_error(
    x = test_x, 
    projection_mat = simpls_opls$projection_mat, 
    xloadings = simpls_opls$X_loadings,
    scale = FALSE,
    Xcenter = simpls_opls$transf$Xcenter,
    Xscale = simpls_opls$transf$Xscale,
    scale_back = TRUE
  )
  
  # Tests
  expect_lt(
    sum(abs(simpls_projection$weights - simpls_basics$weights)), 
    tol,
    label = "SIMPLS scale=FALSE: opls_for_projection vs opls_get_basics weights"
  )
  
  expect_lt(
    sum(abs(simpls_basics$weights - simpls_opls$weights)), 
    tol,
    label = "SIMPLS scale=FALSE: opls_get_basics vs opls weights"
  )
  
  expect_lt(
    sum(abs(simpls_basics_predict - simpls_opls_predict)), 
    tol2,
    label = "SIMPLS scale=FALSE: predict_opls consistency"
  )
  
  expect_lt(
    sum(abs(simpls_projection$scores - simpls_opls$scores)), 
    tol,
    label = "SIMPLS scale=FALSE: opls_for_projection vs opls scores"
  )
  
  expect_lt(
    sum(abs(simpls_gesearch$pred_response - simpls_opls_predict[, ncomp_test])), 
    tol2,
    label = "SIMPLS scale=FALSE: opls_gesearch vs opls predictions"
  )
  
  expect_lt(
    abs(drop(simpls_gesearch$rmse_reconstruction) - drop(simpls_opls_rec)), 
    tol,
    label = "SIMPLS scale=FALSE: opls_gesearch vs reconstruction_error"
  )
  
  # Xscale should be all ones when scale = FALSE
  expect_true(
    all(simpls_basics$transf$Xscale == 1),
    label = "SIMPLS scale=FALSE: Xscale should be all ones (opls_get_basics)"
  )
  
  expect_true(
    all(simpls_opls$transf$Xscale == 1),
    label = "SIMPLS scale=FALSE: Xscale should be all ones (opls)"
  )
})


test_that("mpls implementations produce consistent results with scale = TRUE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40
  tol <- 1e-8
  
  # Fit models
  mpls_basics <- resemble:::opls_get_basics(
    X = train_x, 
    Y = as.matrix(train_y), 
    ncomp = ncomp_test,
    scale = TRUE,
    maxiter = 100,
    tol = 1e-06, 
    algorithm = "mpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  mpls_opls <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "mpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  mpls_projection <- resemble:::opls_for_projection(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,            
    maxiter = 1000,
    tol = 1e-06,
    pcSelmethod = "manual",
    pcSelvalue = ncomp_test, 
    algorithm = "mpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  mpls_gesearch <- resemble:::opls_gesearch(
    Xr = train_x, 
    Yr = as.matrix(train_y), 
    Xu = test_x, 
    ncomp = ncomp_test,
    scale = TRUE,     
    response = TRUE, 
    reconstruction = TRUE,
    similarity = TRUE,
    fresponse = TRUE,
    algorithm = "mpls"
  )
  
  # Predictions
  mpls_basics_predict <- resemble:::predict_opls(
    bo = mpls_basics$bo, 
    b = mpls_basics$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = mpls_basics$transf$Xscale
  )
  
  mpls_opls_predict <- resemble:::predict_opls(
    bo = mpls_opls$bo, 
    b = mpls_opls$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = mpls_opls$transf$Xscale
  )
  
  # Reconstruction error
  mpls_opls_rec <- resemble:::reconstruction_error(
    x = test_x, 
    projection_mat = mpls_opls$projection_mat, 
    xloadings = mpls_opls$X_loadings,
    scale = TRUE,
    Xcenter = mpls_opls$transf$Xcenter,
    Xscale = mpls_opls$transf$Xscale,
    scale_back = TRUE
  )
  
  # Tests
  expect_lt(
    sum(abs(mpls_projection$weights - mpls_basics$weights)), 
    tol,
    label = "mpls scale=TRUE: opls_for_projection vs opls_get_basics weights"
  )
  
  expect_lt(
    sum(abs(mpls_basics$weights - mpls_opls$weights)), 
    tol,
    label = "mpls scale=TRUE: opls_get_basics vs opls weights"
  )
  
  expect_lt(
    sum(abs(mpls_basics_predict - mpls_opls_predict)), 
    tol,
    label = "mpls scale=TRUE: predict_opls consistency"
  )
  
  expect_lt(
    sum(abs(mpls_projection$scores - mpls_opls$scores)), 
    tol,
    label = "mpls scale=TRUE: opls_for_projection vs opls scores"
  )
  
  expect_lt(
    sum(abs(mpls_gesearch$pred_response - mpls_opls_predict[, ncomp_test])), 
    tol,
    label = "mpls scale=TRUE: opls_gesearch vs opls predictions"
  )
  
  expect_lt(
    abs(drop(mpls_gesearch$rmse_reconstruction) - drop(mpls_opls_rec)), 
    tol,
    label = "mpls scale=TRUE: opls_gesearch vs reconstruction_error"
  )
})


test_that("mpls implementations produce consistent results with scale = FALSE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40
  tol <- 1e-8
  
  # Fit models
  mpls_basics <- resemble:::opls_get_basics(
    X = train_x, 
    Y = as.matrix(train_y), 
    ncomp = ncomp_test,
    scale = FALSE,
    maxiter = 100,
    tol = 1e-06, 
    algorithm = "mpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  mpls_opls <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "mpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  mpls_projection <- resemble:::opls_for_projection(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,            
    maxiter = 1000,
    tol = 1e-06,
    pcSelmethod = "manual",
    pcSelvalue = ncomp_test, 
    algorithm = "mpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  mpls_gesearch <- resemble:::opls_gesearch(
    Xr = train_x, 
    Yr = as.matrix(train_y), 
    Xu = test_x, 
    ncomp = ncomp_test,
    scale = FALSE,     
    response = TRUE, 
    reconstruction = TRUE,
    similarity = TRUE,
    fresponse = TRUE,
    algorithm = "mpls"
  )
  
  # Predictions
  mpls_basics_predict <- resemble:::predict_opls(
    bo = mpls_basics$bo, 
    b = mpls_basics$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = mpls_basics$transf$Xscale
  )
  
  mpls_opls_predict <- resemble:::predict_opls(
    bo = mpls_opls$bo, 
    b = mpls_opls$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = mpls_opls$transf$Xscale
  )
  
  # Reconstruction error
  mpls_opls_rec <- resemble:::reconstruction_error(
    x = test_x, 
    projection_mat = mpls_opls$projection_mat, 
    xloadings = mpls_opls$X_loadings,
    scale = FALSE,
    Xcenter = mpls_opls$transf$Xcenter,
    Xscale = mpls_opls$transf$Xscale,
    scale_back = TRUE
  )
  
  # Tests
  expect_lt(
    sum(abs(mpls_projection$weights - mpls_basics$weights)), 
    tol,
    label = "mpls scale=FALSE: opls_for_projection vs opls_get_basics weights"
  )
  
  expect_lt(
    sum(abs(mpls_basics$weights - mpls_opls$weights)), 
    tol,
    label = "mpls scale=FALSE: opls_get_basics vs opls weights"
  )
  
  expect_lt(
    sum(abs(mpls_basics_predict - mpls_opls_predict)), 
    tol,
    label = "mpls scale=FALSE: predict_opls consistency"
  )
  
  expect_lt(
    sum(abs(mpls_projection$scores - mpls_opls$scores)), 
    tol,
    label = "mpls scale=FALSE: opls_for_projection vs opls scores"
  )
  
  expect_lt(
    sum(abs(mpls_gesearch$pred_response - mpls_opls_predict[, ncomp_test])), 
    tol,
    label = "mpls scale=FALSE: opls_gesearch vs opls predictions"
  )
  
  expect_lt(
    abs(drop(mpls_gesearch$rmse_reconstruction) - drop(mpls_opls_rec)), 
    tol,
    label = "mpls scale=FALSE: opls_gesearch vs reconstruction_error"
  )
  
  # Xscale should be all ones when scale = FALSE
  expect_true(
    all(mpls_basics$transf$Xscale == 1),
    label = "mpls scale=FALSE: Xscale should be all ones (opls_get_basics)"
  )
  
  expect_true(
    all(mpls_opls$transf$Xscale == 1),
    label = "mpls scale=FALSE: Xscale should be all ones (opls)"
  )
})


test_that("PLS and SIMPLS produce equivalent predictions with scale = TRUE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40
  tol2 <- tol <- 1e-5  # Looser tolerance for cross-algorithm comparison
  tol3 <- 1e-3  # Looser tolerance for coefficients comparison
  # Fit PLS (NIPALS)
  pls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  # Fit SIMPLS
  simpls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  # Predictions
  pls_predict <- resemble:::predict_opls(
    bo = pls_fit$bo, 
    b = pls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = pls_fit$transf$Xscale
  )
  
  simpls_predict <- resemble:::predict_opls(
    bo = simpls_fit$bo, 
    b = simpls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = TRUE,
    Xscale = simpls_fit$transf$Xscale
  )
  
  # PLS and SIMPLS should produce equivalent predictions
  expect_lt(
    sum(abs(pls_predict - simpls_predict)), 
    tol2,
    label = "PLS vs SIMPLS scale=TRUE: predictions equivalence"
  )
  
  # Coefficients should also match
  expect_lt(
    sum(abs(pls_fit$coefficients - simpls_fit$coefficients)), 
    tol3,
    label = "PLS vs SIMPLS scale=TRUE: coefficients equivalence"
  )
})


test_that("PLS and SIMPLS produce equivalent predictions with scale = FALSE", {
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  tol2 <- 1e-5
  tol3 <- 1e-3 
  # Preprocessing
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det

  # Split data
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  test_idx <- NIRsoil$train == 0 & !is.na(NIRsoil$Ciso)
  
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  test_x <- NIRsoil$spc_pr[test_idx, ]
  
  ncomp_test <- 40

  # Fit PLS (NIPALS)
  pls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  # Fit SIMPLS
  simpls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  # Predictions
  pls_predict <- resemble:::predict_opls(
    bo = pls_fit$bo, 
    b = pls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = pls_fit$transf$Xscale
  )
  
  simpls_predict <- resemble:::predict_opls(
    bo = simpls_fit$bo, 
    b = simpls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = test_x,
    scale = FALSE,
    Xscale = simpls_fit$transf$Xscale
  )
  
  # PLS and SIMPLS should produce equivalent predictions
  expect_lt(
    sum(abs(pls_predict - simpls_predict)), 
    tol2,
    label = "PLS vs SIMPLS scale=FALSE: predictions equivalence"
  )
  
  # Coefficients should also match
  expect_lt(
    sum(abs(pls_fit$coefficients - simpls_fit$coefficients)), 
    tol3,
    label = "PLS vs SIMPLS scale=FALSE: coefficients equivalence"
  )
})


test_that("PLS overfits training data with many components (scale = TRUE)", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 0, p = 1, w = 17
  )
  NIRsoil$spc_pr <- sg_det
  
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  
  ncomp_test <- 200 #min(dim(train_x)) - 1
  
  # Fit PLS (NIPALS)
  pls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_pred <- resemble:::predict_opls(
    bo = pls_fit$bo, 
    b = pls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = train_x,
    scale = TRUE,
    Xscale = pls_fit$transf$Xscale
  )
  
  pred <- pls_pred[, ncomp_test]
  obs <- train_y
  residuals <- pred - obs
  # With 100 components, training predictions should nearly perfectly match observations
  # Any systematic bias or corruption would show up here
  
  # Mean bias should be essentially zero (not shifted up or down)
  expect_lt(abs(mean(residuals)), 1e-5, 
            label = "PLS scale=TRUE: no systematic bias")
  
  # use relative error
  expect_lt(max(abs(residuals)) / diff(range(obs)), 0.05, 
            label = "PLS scale=TRUE: max relative error < 5%")
  
  # R2 should be essentially 1
  expect_gt(cor(pred, obs)^2, 0.98, 
            label = "PLS scale=TRUE: near-perfect fit")
  
  # Slope of pred vs obs should be ~1 (no scaling issues)
  slope <- cov(pred, obs) / var(obs)
  expect_lt(abs(slope - 1), 0.98, 
            label = "PLS scale=TRUE: slope near 1")
  
  # Fit SIMPLS
  simpls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = TRUE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_pred <- resemble:::predict_opls(
    bo = simpls_fit$bo, 
    b = simpls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = train_x,
    scale = TRUE,
    Xscale = simpls_fit$transf$Xscale
  )
  
  pred_simpls <- simpls_pred[, ncomp_test]
  residuals_simpls <- pred_simpls - obs
  
  expect_lt(abs(mean(residuals_simpls)), 1e-5, 
            label = "SIMPLS scale=TRUE: no systematic bias")
  expect_lt(max(abs(residuals_simpls)) / diff(range(obs)), 0.05, 
            label = "SIMPLS scale=TRUE: max relative error < 5%")
  expect_gt(cor(pred_simpls, obs)^2, 0.98, 
            label = "SIMPLS scale=TRUE: near-perfect fit")
  
  slope_simpls <- cov(pred_simpls, obs) / var(obs)
  expect_lt(abs(slope_simpls - 1), 0.98, 
            label = "SIMPLS scale=TRUE: slope near 1")
  
  # PLS and SIMPLS should be equivalent
  expect_lt(max(abs(pred - pred_simpls)), 1e-5,
            label = "PLS vs SIMPLS scale=TRUE: equivalent predictions")
})


test_that("PLS overfits training data with many components (scale = FALSE)", {
  skip_on_cran()
  skip_if_not_installed("prospectr")
  
  library(prospectr)
  data(NIRsoil, package = "prospectr")
  
  sg_det <- savitzkyGolay(
    detrend(NIRsoil$spc, wav = as.numeric(colnames(NIRsoil$spc))),
    m = 1, p = 1, w = 7
  )
  NIRsoil$spc_pr <- sg_det
  
  train_idx <- NIRsoil$train == 1 & !is.na(NIRsoil$Ciso)
  train_x <- NIRsoil$spc_pr[train_idx, ]
  train_y <- NIRsoil$Ciso[train_idx]
  
  ncomp_test <- floor(0.75 * min(dim(train_x)) - 1)
  
  # Fit PLS (NIPALS)
  pls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "pls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  pls_pred <- resemble:::predict_opls(
    bo = pls_fit$bo, 
    b = pls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = train_x,
    scale = FALSE,
    Xscale = pls_fit$transf$Xscale
  )
  
  pred <- pls_pred[, ncomp_test]
  obs <- train_y
  residuals <- pred - obs
  # plot(pred, obs)
  # abline(0, 1)
  
  expect_lt(abs(mean(residuals)), 1e-5, 
            label = "PLS scale=FALSE: no systematic bias")
  expect_lt(max(abs(residuals)) / diff(range(obs)), 0.05, 
            label = "PLS scale=FALSE: max relative error < 5%")
  expect_gt(cor(pred, obs)^2, 0.9999, 
            label = "PLS scale=FALSE: near-perfect fit")
  
  slope <- cov(pred, obs) / var(obs)
  expect_lt(abs(slope - 1), 1e-4, 
            label = "PLS scale=FALSE: slope near 1")
  
  # Fit SIMPLS
  simpls_fit <- resemble:::opls(
    X = train_x, 
    Y = as.matrix(train_y),
    ncomp = ncomp_test,
    scale = FALSE,           
    maxiter = 1000,
    tol = 1e-06,
    algorithm = "simpls", 
    xls_min_w = 3, 
    xls_max_w = 15
  )
  
  simpls_pred <- resemble:::predict_opls(
    bo = simpls_fit$bo, 
    b = simpls_fit$coefficients, 
    ncomp = ncomp_test, 
    newdata = train_x,
    scale = FALSE,
    Xscale = simpls_fit$transf$Xscale
  )
  
  pred_simpls <- simpls_pred[, ncomp_test]
  residuals_simpls <- pred_simpls - obs
  
  expect_lt(abs(mean(residuals_simpls)), 1e-5, 
            label = "SIMPLS scale=FALSE: no systematic bias")
  expect_lt(max(abs(residuals)) / diff(range(obs)), 0.05, 
            label = "SIMPLS scale=FALSE: max relative error < 5%")
  expect_gt(cor(pred_simpls, obs)^2, 0.9999, 
            label = "SIMPLS scale=FALSE: near-perfect fit")
  
  slope_simpls <- cov(pred_simpls, obs) / var(obs)
  expect_lt(abs(slope_simpls - 1), 1e-3, 
            label = "SIMPLS scale=FALSE: slope near 1")
  
  # PLS and SIMPLS should be equivalent
  expect_lt(mean(abs(pred - pred_simpls)), 0.04,
            label = "PLS vs SIMPLS scale=FALSE: equivalent predictions")
})
