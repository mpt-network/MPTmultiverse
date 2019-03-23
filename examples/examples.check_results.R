load(file = system.file("extdata", "results_bayen_kuhlmann.RData",
                        package = "MPTmultiverse"))
## prints checks to console
check_results(results)

## returns tibble with single row
check_set(results)
