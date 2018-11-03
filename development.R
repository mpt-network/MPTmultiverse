
library("devtools")
library("testthat")
load_all()
test()

document()

check()

library(usethis) ## see: https://github.com/r-lib/usethis

use_test("mptinr")
