
# ------------------------------------------------------------------------------
# MPT model definition & Data

EQN_FILE <- system.file("extdata", "prospective_memory.eqn", package = "MPTmultiverse")
DATA_FILE <- system.file("extdata", "smith_et_al_2011.csv", package = "MPTmultiverse")

### if .csv format uses semicolons ";" (e.g., German format):
# data <- read.csv2(DATA_FILE, fileEncoding = "UTF-8-BOM")
### if .csv format uses commata "," (international format):
data <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM")
data <- data[c(1:10, 113:122),]
head(data)

COL_CONDITION <- "WM_EX"  # name of the variable encoding group membership

# experimental condition should be labeled meaningfully ----
unique(data[[COL_CONDITION]])

data[[COL_CONDITION]] <- factor(
  data[[COL_CONDITION]]
  , levels = 1:2
  , labels = c("low_WM", "high_WM")
)

# set test options for a quick and unreliable run:
mpt_options("test")
mpt_options() # to

all_supported_methods <- c(
  "asymptotic_complete"
  , "asymptotic_no"
  , "pb_no"
  , "npb_no"
  , "simple"
  , "simple_pooling"
  , "trait"
  , "trait_uncorrelated"
  , "beta"
)


only_asymptotic <- fit_mpt(
  method = "asymptotic_no"
  , dataset = DATA_FILE
  , data = data
  , model = EQN_FILE
  , condition = COL_CONDITION
)

dplyr::glimpse(only_asymptotic)

all_bootstrap <- fit_mpt(
  method = c("pb_no", "npb_no")
  , dataset = DATA_FILE
  , data = data
  , model = EQN_FILE
  , condition = COL_CONDITION
)

dplyr::glimpse(all_bootstrap)
