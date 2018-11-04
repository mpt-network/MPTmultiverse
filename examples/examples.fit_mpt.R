
# ------------------------------------------------------------------------------
# MPT model definition & Data

EQN_FILE <- system.file("extdata", "prospective_memory.eqn", package = "MPTmultiverse")
DATA_FILE <- system.file("extdata", "smith_et_al_2011.csv", package = "MPTmultiverse")

### if .csv format uses semicolons ";" (e.g., German format):
# data <- read.csv2(DATA_FILE, fileEncoding = "UTF-8-BOM")
### if .csv format uses commata "," (international format):
data <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM")
data <- data[c(1:10, 113:122),]  ## select only subset of data for example
head(data)

COL_CONDITION <- "WM_EX"  # name of the variable encoding group membership

# experimental condition should be labeled meaningfully ----
unique(data[[COL_CONDITION]])

data[[COL_CONDITION]] <- factor(
  data[[COL_CONDITION]]
  , levels = 1:2
  , labels = c("low_WM", "high_WM")
)

# define core parameters:
CORE <- c("C1", "C2")

## save options so they can be reset later:
op <- mpt_options() 

# set test options for a quick and unreliable run:
mpt_options("test")
mpt_options("n.CPU" = 1) # use 1 core to run on CRAN

\dontrun{
## to reset default options (which you would want) use:
mpt_options("default")

mpt_options() # to see the settings 
## Note: settings are also saved in the results tibble
  
## without specifying method, all are used per default
fit_all <- fit_mpt(
  dataset = DATA_FILE
  , data = data
  , model = EQN_FILE
  , condition = COL_CONDITION
  , core = CORE
)

### Analysis of results requires dplyr and tidyr (or just 'tidyverse')
library("dplyr")
library("tidyr")

glimpse(fit_all) 
## first few columns identify model, data, and estimation approach/method
## remaining columns are list columns containing the results for each method
## use unnest to work with each of the results columns

fit_all %>% 
  select(method, pooling, est_group) %>% 
  unnest() 


fit_all %>% 
  select(method, pooling, gof_group) %>% 
  unnest() %>% 
  as.data.frame()

  
}

### Also possible to only use individual methods:

only_asymptotic <- fit_mpt(
  method = "asymptotic_no"
  , dataset = DATA_FILE
  , data = data
  , model = EQN_FILE
  , condition = COL_CONDITION
  , core = CORE
)

dplyr::glimpse(only_asymptotic)

all_bootstrap <- fit_mpt(
  method = c("pb_no", "npb_no")
  , dataset = DATA_FILE
  , data = data
  , model = EQN_FILE
  , condition = COL_CONDITION
  , core = CORE
)

dplyr::glimpse(all_bootstrap)

## reset default options:
mpt_options("default")
