#' Options Settings for MPT Comparison
#' 
#' Set and examine a variety of \emph{options} which affect the way MPT models
#' are estimated.
#' 
#' @examples 
#' # Examine options:
#' mpt_options()
#' 
#' # Set number of MCMC chains to 20:
#' mpt_options(n.chain = 20)
#' mpt_options()
#' 
#' @export
 
mpt_options <- function(...){
  
  fetched <- getOption("mpt.comparison")
  args <- c(...)
  
  if(length(args)==0L) return(fetched)
  if(args=="test"){
    changed <- set_test()
  } else {
    changed <- lapply(
      X = fetched
      , FUN = function(x, args){
        sub_args <- args[names(args)%in%names(x)]
        x[names(sub_args)] <- sub_args
        x
      }
      , args = args
    )
    sub_args <- args[names(args)%in%names(fetched)]
    changed[names(sub_args)] <- sub_args
  }
  options(list(mpt.comparison = changed))
}



#' @keywords internal

set_test <- function() { # nocov start
  cat("Setting options for a quick test run.\nDo not interpret results!")
  
  list(
    mptinr = list(
      bootstrap_samples = 2e1
      , n.optim = 2
      , n.CPU = parallel::detectCores()
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
      , n.CPU = parallel::detectCores()
    )
    , ci_size = c(.025, .1, .9, .975)
    , max_ci_indiv = .99
    , save_models = FALSE
  )
}