
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

\dontrun{
op <- mpt_options() 
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

mpt_options(op) ## reset options  
}

load(system.file("extdata", "prospective_memory_example.rda", package = "MPTmultiverse"))

# Although we requested all 10 methods, only 9 worked:
fit_all$method
# Jags variant of beta MPT is missing.

# the returned method has a plot method. For example, for the group-level estimates:
plot(fit_all, which = "est")

\dontrun{
### Full analysis of results requires dplyr and tidyr (or just 'tidyverse')
library("dplyr")
library("tidyr")

## first few columns identify model, data, and estimation approach/method
## remaining columns are list columns containing the results for each method
## use unnest to work with each of the results columns
glimpse(fit_all) 

## Let us inspect the group-level estimates
fit_all %>% 
  select(method, pooling, est_group) %>% 
  unnest() 

## which we can plot again
plot(fit_all, which = "est")

## Next we take a look at the GoF
fit_all %>% 
  select(method, pooling, gof_group) %>% 
  unnest() %>% 
  as.data.frame()

# Again, we can plot it as well
plot(fit_all, which = "gof2")  ## use "gof1" for overall GoF

## Finally, we take a look at the differences between conditions
fit_all %>% 
  select(method, pooling, test_between) %>% 
  unnest() 

# and then we plot it
plot(fit_all, which = "test_between")


### Also possible to only use individual methods:
only_asymptotic <- fit_mpt(
  method = "asymptotic_no"
  , dataset = DATA_FILE
  , data = data
  , model = EQN_FILE
  , condition = COL_CONDITION
  , core = CORE
)

glimpse(only_asymptotic)

bayes_complete <- fit_mpt(
  method = c("simple_pooling")
  , dataset = DATA_FILE
  , data = data
  , model = EQN_FILE
  , condition = COL_CONDITION
  , core = CORE
)
glimpse(bayes_complete)

}