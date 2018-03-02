#' @export

set_meticulous <- function() { # nocov start
  
  op <- options()
  op_mpt <- list(
    mpt.comparison = list(
      mptinr = list(
        bootstrap_samples = 2e3
        , n.optim = 1e2
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
    )
  )
  
  toset <- !(names(op_mpt) %in% names(op))
  if(any(toset)) options(op_mpt[toset])
  
  
  invisible()
}

#' @export

set_test <- function() { # nocov start
  
  op <- options()
  op_mpt <- list(
    mpt.comparison = list(
      mptinr = list(
        bootstrap_samples = 2e1
        , n.optim = 5
        , n.CPU = NULL
      )
      , treebugs = list(
        n.chain = 4
        , n.iter = 2e2
        , n.adapt = 1e1
        , n.burnin = 1e2
        , n.thin = 1e0
        , Rhat_max = 10
        , Neff_min = 1
        , extend_max = 0
        , n.PPP = 20
        , n.CPU = 4
      )
      , ci_size = c(.025, .1, .9, .975)
      , max_ci_indiv = .99
    )
  )
  
  options(op_mpt)
  
  
  invisible()
}