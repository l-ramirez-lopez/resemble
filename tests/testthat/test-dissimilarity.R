context("test-dissimilarity")



test_that("dissimilarity works", {
  nirdata <- data("NIRsoil", package = "prospectr")
  
  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
  
  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]
  
  Xu <- Xu[!is.na(Yu), ][1:20, ]
  Xr <- Xr[!is.na(Yr), ][1:40, ]
  
  Yu <- Yu[!is.na(Yu)][1:20]
  Yr <- Yr[!is.na(Yr)][1:40]
  
  dsm_pca <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = c("pca"),
    Yr = Yr, gh = TRUE, pc_selection = list("opc", 15),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_comp <- 6
  dsm_pls <- dissimilarity(
    Xr = standardNormalVariate(Xr), Xu = standardNormalVariate(Xu),
    diss_method = c("pls"),
    Yr = Yr, gh = TRUE, pc_selection = list("opc", 15),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_pls <- 15

  dsm_pca_var <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = c("pca"),
    Yr = Yr, gh = TRUE, pc_selection = list("var", 0.01),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_comp_var <- 2
  dsm_pls_var <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = c("pls"),
    Yr = Yr, gh = TRUE, pc_selection = list("var", 0.01),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_pls_var <- 2
  
  
  dsm_euclid <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = "euclid",
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  
  dsm_euclid_xu <- dissimilarity(
    Xr = Xu[1:10, ],
    diss_method = "euclid",
    center = FALSE, scale = FALSE
  )$dissimilarity
  dsm_euclid_xu <- ((dsm_euclid_xu^2) * ncol(Xu))^0.5
  dist_euclid_xu <- as.matrix(dist(Xu[1:10, ]))
  
  dsm_cor <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = "cor",
    return_projection = TRUE,
    ws = 11,
    center = TRUE, scale = FALSE
  )
  
  output_names_pca <- names(dsm_pca)
  output_names_pls <- names(dsm_pls)
  expected_names <- c("dissimilarity", "projection", "gh", "documentation")
  
  expect_is(dsm_pca, "list")
  expect_is(dsm_pls, "list")
  expect_is(dsm_pca_var, "list")
  expect_is(dsm_pls_var, "list")
  expect_is(dsm_euclid, "list")
  expect_is(dsm_euclid, "list")
  expect_true(dsm_pca$projection$n_components == expected_n_comp)
  expect_true(dsm_pls$projection$n_components == expected_n_pls)
  expect_true(dsm_pca_var$projection$n_components == expected_n_comp_var)
  expect_true(dsm_pls_var$projection$n_components == expected_n_pls_var)
  expect_true(all(expected_names %in% output_names_pca))
  expect_true(all(expected_names %in% output_names_pls))
  expect_true(sum(abs(round(dist_euclid_xu - dsm_euclid_xu, 5))) == 0)
})














test_that("dissimilarity large sets works", {
  skip_on_cran()
  skip_on_travis()
  nirdata <- data("NIRsoil", package = "prospectr")

  Xu <- NIRsoil$spc[!as.logical(NIRsoil$train), ]
  Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]

  Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
  Xr <- NIRsoil$spc[as.logical(NIRsoil$train), ]

  Xu <- Xu[!is.na(Yu), ]
  Xr <- Xr[!is.na(Yr), ]

  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]

  dsm_pca <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = c("pca"),
    Yr = Yr, gh = TRUE, pc_selection = list("opc", 30),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_comp <- 24
  dsm_pls <- dissimilarity(
    Xr = standardNormalVariate(Xr), Xu = standardNormalVariate(Xu),
    diss_method = c("pls"),
    Yr = Yr, gh = TRUE, pc_selection = list("opc", 30),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_pls <- 10

  dsm_pca_var <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = c("pca"),
    Yr = Yr, gh = TRUE, pc_selection = list("var", 0.02),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_comp_var <- 2
  dsm_pls_var <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = c("pls"),
    Yr = Yr, gh = TRUE, pc_selection = list("var", 0.02),
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )
  expected_n_pls_var <- 2


  dsm_euclid <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = "euclid",
    return_projection = TRUE,
    center = TRUE, scale = TRUE
  )

  dsm_euclid_xu <- dissimilarity(
    Xr = Xu[1:10, ],
    diss_method = "euclid",
    center = FALSE, scale = FALSE
  )$dissimilarity
  dsm_euclid_xu <- ((dsm_euclid_xu^2) * ncol(Xu))^0.5
  dist_euclid_xu <- as.matrix(dist(Xu[1:10, ]))

  dsm_cor <- dissimilarity(
    Xr = Xr, Xu = Xu,
    diss_method = "cor",
    return_projection = TRUE,
    ws = 11,
    center = TRUE, scale = FALSE
  )

  output_names_pca <- names(dsm_pca)
  output_names_pls <- names(dsm_pls)
  expected_names <- c("dissimilarity", "projection", "gh", "documentation")

  expect_is(dsm_pca, "list")
  expect_is(dsm_pls, "list")
  expect_is(dsm_pca_var, "list")
  expect_is(dsm_pls_var, "list")
  expect_is(dsm_euclid, "list")
  expect_is(dsm_euclid, "list")
  expect_true(dsm_pca$projection$n_components == expected_n_comp)
  expect_true(dsm_pls$projection$n_components == expected_n_pls)
  expect_true(dsm_pca_var$projection$n_components == expected_n_comp_var)
  expect_true(dsm_pls_var$projection$n_components == expected_n_pls_var)
  expect_true(all(expected_names %in% output_names_pca))
  expect_true(all(expected_names %in% output_names_pls))
  expect_true(sum(abs(round(dist_euclid_xu - dsm_euclid_xu, 5))) == 0)
})
