context("MPTinR basic tests")

test_that("No-pooling approaches work", {
  
  testthat::skip_on_cran()
  
  EQN_FILE <- system.file("extdata", "prospective_memory.eqn", package = "MPTmultiverse")
  DATA_FILE <- system.file("extdata", "smith_et_al_2011.csv", package = "MPTmultiverse")
  
  data <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM")
  data <- data[c(1:5, 113:118),]
  COL_CONDITION <- "WM_EX"
  data[[COL_CONDITION]] <- factor(
    data[[COL_CONDITION]]
    , levels = 1:2
    , labels = c("low_WM", "high_WM")
  )
  op <- mpt_options()
  capture_output(mpt_options("test"))
  
  only_asymptotic <- fit_mpt(
    method = "asymptotic_no"
    , dataset = DATA_FILE
    , data = data
    , model = EQN_FILE
    , condition = COL_CONDITION
  )
  expect_equal(nrow(only_asymptotic), 1)
  expect_equal(only_asymptotic$pooling, "no")
  expect_equal(only_asymptotic$method, "asymptotic")
  
  only_pb <- fit_mpt(
    method = "pb_no"
    , dataset = DATA_FILE
    , data = data
    , model = EQN_FILE
    , condition = COL_CONDITION
  )
  expect_equal(nrow(only_pb), 1)
  expect_equal(only_pb$pooling, "no")
  expect_equal(only_pb$method, "PB/MLE")
  
  only_npb <- fit_mpt(
    method = "npb_no"
    , dataset = DATA_FILE
    , data = data
    , model = EQN_FILE
    , condition = COL_CONDITION
  )
  expect_equal(nrow(only_npb), 1)
  expect_equal(only_npb$pooling, "no")
  expect_equal(only_npb$method, "NPB/MLE")
  
})
