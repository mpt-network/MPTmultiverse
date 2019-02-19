
library("devtools")
library("testthat")
load_all()

test()

document()

check()

build_vignettes()
build()

library(usethis) ## see: https://github.com/r-lib/usethis

use_test("mptinr")

use_build_ignore("development.R")
use_package("limSolve", "Imports")

#### stuff for main example below

save(fit_all, file = "inst/extdata/prospective_memory_example.rda")
