#' Main function
#'
#' Does a lot of nice stuff
#'
#' @param method Character. A vector specifying which analysis approaches should be performed (see Description below).
#'   Defaults to all available methods.
#' @param dataset Character. The name of the dataset to be analyzed. 
#' @param data A \code{data.frame} containing the data.
#' @param model A model definition, typically the path to an \code{.eqn} file.
#' @param id Character. Name of the column that contains the subject identifier.
#'   If not specified, it is assumed that each row represents observations from one participant.
#' @param condition Character. Name of the column specifying a between-subjects factor.
#'   If not specified, no between-subjects comparisons are performed.
#' @param core character vector defining the core parameters of interest, e.g.,
#'   \code{core = c("Dn", "Do")}. All other parameters are treated as auxiliary parameters.
#' @examples examples/examples.fit_mpt.R
#' 
#' @details 
#' Maximum-likelihood estimation with MPTinR:
#' \itemize{
#'   \item{\code{"asymptotic_complete"}: }{Asymptotic ML theory, complete pooling}
#'   \item{\code{"asymptotic_no"}: }{ Asymptotic ML theory, no pooling}
#'   \item{\code{"pb_no"}: }{Parametric bootstrap, no pooling}
#'   \item{\code{"npb_no"}: }{Nonparametric bootstrap, no pooling}
#' }
#'  
#' Bayesian estimation with TreeBUGS
#' \itemize{
#'   \item{\code{"simple"}: }{Bayesian estimation, no pooling (C++, \link[TreeBUGS]{simpleMPT})}
#'   \item{\code{"simple_pooling"}: }{Bayesian estimation, complete pooling (C++, \link[TreeBUGS]{simpleMPT})}
#'   \item{\code{"trait"}: }{latent-trait model, partial pooling (JAGS, \link[TreeBUGS]{traitMPT})}
#'   \item{\code{"trait_uncorrelated"}: }{latent-trait model without correlation parameters, partial pooling (JAGS, \link[TreeBUGS]{traitMPT})}
#'   \item{\code{"beta"}: }{beta-MPT model, partial pooling (JAGS, \link[TreeBUGS]{betaMPT})}
#'   \item{\code{"betacpp"}: }{beta-MPT model, partial pooling (C++, \link[TreeBUGS]{betaMPTcpp})}
#' }
#'
#' @export

fit_mpt <- function(
  method
  , dataset
  , data
  , model
  , id = NULL
  , condition = NULL
  , core = NULL
) {
  
  available_methods <- c(
    # MPTinR ----
    "asymptotic_complete"
    , "asymptotic_no"
    , "pb_no"
    , "npb_no"
    # TreeBUGS ----
    , "simple"
    , "simple_pooling"
    , "trait"
    , "trait_uncorrelated"
    , "beta"
    , "betacpp"
  )
  
  if(missing(method)) {
    method <- available_methods
  }
  
  method <- match.arg(
    arg = method
    , choices = available_methods
    , several.ok = TRUE
  )
  
  # set options ----
  silent_jags <- getOption("MPTmultiverse")$silent_jags
  runjags::runjags.options(silent.jags = silent_jags, 
                           silent.runjags = silent_jags)
  
  # prepare data ----
  if (missing(data)) {
    data <- as.data.frame(readr::read_csv(dataset))
  }
  
  if(is.null(condition)) {
    data$ExpCond <- "no_condition"
    condition <- "ExpCond"
  }
  
  if(is.null(id)) {
    data$Subject <- 1:nrow(data)
    id <- "Subject"
  }
  
  # Ensure that all variables are character
  data$ExpCond <- as.character(data[[condition]])
  data$Subject <- as.character(data[[id]])
  
  
  # check MPT file
  mpt_model <- TreeBUGS::readEQN(model)
  
  if(!is.data.frame(mpt_model)) {
    "I can't comprehend your .eqn file."
  }

  
  # remove extraneous colums and check if all specified columns are present
  # in data
  freq_cols <- get_eqn_categories(model)
  valid_cols <- c(id, condition, freq_cols)
  check_cols <- valid_cols %in% colnames(data)
  
  if(!all(check_cols)) {
    stop("Variable \"", paste(valid_cols[!check_cols], collapse = ", "), "\" not found in data.frame.")
  }
  
  data <- data[, valid_cols]
  
  
  # Check NAs ----
  nas_found <- unlist(lapply(X = data, FUN = anyNA))
  
  if(any(nas_found)) {
    stop("Variable \"", paste(valid_cols[nas_found], collapse = ", "), "\" contains missing values.")
  }
  
  # Check whether freqencies are integer ----
  not_integer <- unlist(lapply(X = data[, freq_cols], FUN = function(x) {
      any(as.integer(x)!=x)
    }
  ))
  
  if(any(not_integer)) {
    stop("Variable \"", paste(freq_cols[not_integer], collapse = ", "), "\" contains non-integer values.")
  }

  # Sanity check for TreeBUGS options
  if(any(method %in% c("simple", "simple_pooling", "trait", "trait_uncorrelated", "beta", "betacpp"))) {
    opt <- mpt_options()
    n.samples <- (opt$treebugs$n.iter - opt$treebugs$n.burnin) * opt$n.CPU
    if((opt$treebugs$n.iter - opt$treebugs$n.burnin)<=0) {
      stop("Check your mpt_options(): You specified less iterations (n.iter) than burn-in samples (n.burnin).")
    }
    if(n.samples <= 0 | n.samples < opt$treebugs$Neff_min) {
      warning("With your current mpt_options(), it is not possible to obtain the specified number of effective samples.")
    }
  }
  
  
  # Ensure that id and condition are character, also drops unused levels
  data[[id]] <- as.character(data[[id]])
  data[[condition]] <- as.character(data[[condition]])

  
  res <- list()
  
  # MPTinR part ----
  res[["mptinr"]] <- mpt_mptinr(
    dataset = dataset
    , data = data
    , model = model
    , method = intersect(method, c("asymptotic_complete", "asymptotic_no", "pb_no", "npb_no"))
    , id = id
    , condition = condition
    , core = core
  )
  
  # TreeBUGS part ----
  res[["treebugs"]] <- dplyr::bind_rows(
    purrr::map(
      intersect(method, c("simple", "simple_pooling", "trait", "trait_uncorrelated", "beta", "betacpp"))
      , mpt_treebugs_safe
      , dataset = dataset
      , data = data
      , model = model
      , id = id
      , condition = condition
      , core = core
    )
  )

  y <- dplyr::bind_rows(res)
  class(y) <- c("multiverseMPT", class(y))
  y
}