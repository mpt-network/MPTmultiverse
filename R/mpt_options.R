#' Options Settings for MPT Comparison
#' 
#' Set and examine a variety of \emph{options} which affect the way MPT models
#' are estimated.
#' 
#' @param ... Named parameters to set. Possible values are:
#' \describe{
#'   \item{\strong{bootstrap_samples}}{Numeric. The number of bootstrap samples to be drawn for the calculation parametric bootstrap confidence intervals.}
#'   \item{\strong{n.optim}}{Numeric. The number of optimization runs for the models estimated with maximum-likelihood methods.}
#'   \item{\strong{n.chain}}{Numeric. The number of MCMC chains for the Bayesian models.}
#'   \item{\strong{n.adapt}}{Numeric. The number of iterations for adaptation.}
#'   \item{\strong{n.burnin}}{Numeric. The number of burn-in/warm-up iterations.}
#'   \item{\strong{n.iter}}{Numeric. The total number of iterations to be drawn \emph{after} adaptation (including burnin).}
#'   \item{\strong{n.thin}}{Numeric. Thinning interval.}
#'   \item{\strong{Rhat_max}}{Numeric. The maximum rhat.}
#'   \item{\strong{Neff_min}}{Numeric. The minimum number of effective samples you are willing to accept.}
#'   \item{\strong{extend_max}}{Numeric.}
#'   \item{\strong{n.PPP}}{Numeric. The number of posterior predictive samples drawn for the calculation of fit statistics T_1 and T_2.}
#'   \item{\strong{n.CPU}}{Numeric. The number of CPU cores to use. Defaults to the number of available cores on your machine.}
#'   \item{\strong{ci_size}}{Numeric.}
#'   \item{\strong{max_ci_indiv}}{Numeric.}
#'   \item{\strong{silent}}{Logical.}
#'   \item{\strong{save_models}}{Logical.}
#' }
#'   
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
  
  fetched <- getOption("MPTmultiverse")
  args <- c(...)
  
  if(length(args)==0L) return(fetched)
  
  # Provide some shorthand terms:
  if(is.character(args[[1]])){
    changed <- switch(
      args[[1]]
      , test = set_test_options()
      , default = set_default_options()
    )
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
  options(list(MPTmultiverse = changed))
}



#' @keywords internal

set_test_options <- function() { # nocov start
  cat("Setting options for a quick test run.\nDo not interpret results!")
  
  list(
    mptinr = list(
      bootstrap_samples = 4e1
      , n.optim = 2
      , n.CPU = parallel::detectCores()
    )
    , hmmtree = list(
      n.optim = 2
      , max_classes = 3
      , fisher_information = "expected"
    )
    , treebugs = list(
      n.chain = 4
      , n.iter = 8e2
      , n.adapt = 1e2
      , n.burnin = 1e2
      , n.thin = 1e0
      , Rhat_max = 10
      , Neff_min = 2
      , extend_max = 1
      , n.PPP = 4e1
      , n.CPU = parallel::detectCores()
    )
    , silent = FALSE
    , ci_size = c(.025, .1, .9, .975)
    , max_ci_indiv = .99
    , save_models = FALSE
  )
}


#' @keywords internal

set_default_options <- function() {
  
  list(
    mptinr = list(
      bootstrap_samples = 2e3
      , n.optim = 2e1
      , n.CPU = parallel::detectCores()
    )
    , hmmtree = list(
      n.optim = 2e1
      , max_classes = 2e1
      , fisher_information = "expected"
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
    , silent = TRUE
    , ci_size = c(.025, .1, .9, .975)
    , max_ci_indiv = .99
    , save_models = TRUE
  )
}