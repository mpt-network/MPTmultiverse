.onLoad <- function(libname, pkgname) { # nocov start

  op <- options()
  op_mpt <- list(
    mpt.comparison = list(
      mptinr = list(
        bootstrap_samples = 2e3
        , n.optim = 2e1
        , n.CPU = parallel::detectCores()
      )
      , treebugs = list(
        n.chain = 4
        , n.iter = 5e4
        , n.adapt = 1e4
        , n.burnin = 2e4
        , n.thin = 1e1
        , Rhat_max = 1.01
        , Neff_min = 1e3
        , extend_max = 2e1
        , n.PPP = 5e3
        , n.CPU = parallel::detectCores()
      )
      , ci_size = c(.025, .1, .9, .975)
      , max_ci_indiv = .99
      , save_models = TRUE
    )
  )
  
  toset <- !(names(op_mpt) %in% names(op))
  if(any(toset)) options(op_mpt[toset])


  invisible()
} # nocov end
