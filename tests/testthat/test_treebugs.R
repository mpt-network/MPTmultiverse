context("TreeBUGS basic tests")

test_that("Partial Pooling approaches work", {

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
  capture_output(mpt_options("default"))
  mpt_options(n.chains = 2)  ## use 2 chains, hopefully it still runs on CRAN
  mpt_options(Neff_min = 100)
  mpt_options(n.iter = 50000)
  mpt_options(save_models = FALSE)
  
  set.seed(10)  ## for reproducibility
  
  expect_warning(capture_output(res_bayes <- fit_mpt(
    method = "trait", 
    , dataset = DATA_FILE
    , data = data
    , model = EQN_FILE
    , condition = COL_CONDITION
  )), "The adaptation phase of one or more models was not completed in 10000 iterations",
  fixed = TRUE)
  
})
