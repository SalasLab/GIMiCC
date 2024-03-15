test_that("Wrong tumor name", {
  data_env <- new.env(parent = emptyenv())
  data("Capper_test_betas", envir = data_env, package = "GIMiCC")
  Capper_test_betas <- data_env[["Capper_test_betas"]]
  expect_error(GIMiCC_Deconvo(Capper_test_betas, tumor.type = "glioblastoma"))
  expect_error(GIMiCC_Deconvo(Capper_test_betas, tumor.type = 6))
})

test_that("Wrong input for h", {
  data_env <- new.env(parent = emptyenv())
  data("Capper_test_betas", envir = data_env, package = "GIMiCC")
  Capper_test_betas <- data_env[["Capper_test_betas"]]
  expect_error(GIMiCC_Deconvo(Capper_test_betas, tumor.type = "GBM", h = 10))
  expect_error(GIMiCC_Deconvo(Capper_test_betas, tumor.type = "GBM", h = "hello"))
})

test_that("Not using beta values", {
  data_env <- new.env(parent = emptyenv())
  data("Capper_test_betas", envir = data_env, package = "GIMiCC")
  Capper_test_betas <- data_env[["Capper_test_betas"]]
  Capper_test_mvals <- minfi::logit2(Capper_test_betas)
  expect_error(GIMiCC_Deconvo(Capper_test_mvals, tumor.type = "GBM", h = 5))
})
