#' Multiverse Analysis for MPT Models
#'
#' Performs a multiverse analysis for multinomial processing tree (MPT) models
#' across maximum-likelihood/frequentist and Bayesian estimation approaches. For
#' the frequentist approaches, no pooling (with and without parametric or
#' nonparametric bootstrap) and complete pooling  are implemented using
#' \pkg{MPTinR}. For the Bayesian approaches, no pooling, complete pooling, and
#' three different variants of partial pooling are implemented using
#' \pkg{TreeBUGS}. Requires \code{data} on a by-participant level with each row
#' corresponding to data from one participant (i.e., different response
#' categories correspond to different columns) and the data can contain a single
#' between-subjects condition. Model equations need to be passed as a
#' \code{.eqn} model file and category labels (first column in \code{.eqn} file)
#' need to match the column names in \code{data}. Results are returned in one
#' \code{tibble} with one row per estimation method.
#'
#' @param method \code{character} vector specifying which analysis approaches
#'   should be performed (see Description below). Defaults to all available
#'   methods.
#' @param dataset scalar \code{character} vector. Name of the data set that will
#'   be copied to the results \code{tibble}.
#' @param data A \code{data.frame} containing the data. Column
#'   names need to match category names in \code{model} (i.e., different from
#'   \pkg{MPTinR} behavior, order of categories is not important, matching is
#'   done via name).
#' @param model A model definition, typically the path to an \code{.eqn} model
#'   file containing the model equations. Category names need to match column
#'   names in \code{data}.
#' @param id scalar \code{character} vector. Name of the column that contains
#'   the subject identifier. If not specified, it is assumed that each row
#'   represents observations from one participant.
#' @param condition scalar \code{character} vector. Name of the column
#'   specifying a between-subjects factor. If not specified, no between-subjects
#'   comparisons are performed.
#' @param core \code{character} vector defining the core parameters of interest,
#'   e.g., \code{core = c("Dn", "Do")}. All other parameters are treated as
#'   auxiliary parameters.
#' @example examples/examples.fit_mpt.R
#' 
#' @details This functions is a fancy wrapper for packages \pkg{MPTinR} and
#'   \pkg{TreeBUGS} applying various frequentist and Bayesian estimation methods
#'   to the same data set using a single MPT model and collecting the results
#'   in one \code{tibble} where each row corresponds to one
#'   estimation method. Note that parameter restrictions (e.g., equating
#'   different parameters or fixing them to a constant) need to be part of the
#'   model (i.e., the \code{.eqn} file) and cannot be passed as an argument.
#'   
#'   The settings for the various methods are specified via function
#'   \code{\link{mpt_options}}. The default settings use all available cores for
#'   calculating the boostrap distribution as well as independent MCMC chains
#'   and should be approproaite for most situations.
#'
#'   The data can have a single between-subjects condition (specified via
#'   \code{condition}). This condition can have more than two levels. If
#'   specified, the pairwise differences between each level, the standard error
#'   of the differences, and confidence-intervals of the differences are
#'   calculated for each parameter. Please note that \code{condition} is
#'   silently converted to \code{character} in the output. Thus, a specific
#'   ordering of the \code{factor} levels in the output cannot be guaranteed.
#'   
#'   Parameter differences or other support for within-subject conditions is not
#'   provided. The best course of action for within-subjects conditions is to
#'   simply include separate trees and separate sets of parameters for each
#'   within-subjects condition. This allows to at least compare the estimates
#'   for each within-subjects condition across estimation method.
#'   
#'   \subsection{Implemented Methods}{
#'     Maximum-likelihood estimation with \pkg{MPTinR} via
#'     \code{\link[MPTinR]{fit.mpt}}:
#'     \itemize{
#'       \item{\code{"asymptotic_complete"}: }{Asymptotic ML theory, complete
#'       pooling}
#'       \item{\code{"asymptotic_no"}: }{ Asymptotic ML theory, no pooling}
#'       \item{\code{"pb_no"}: }{Parametric bootstrap, no pooling}
#'       \item{\code{"npb_no"}: }{Nonparametric bootstrap, no pooling}
#'     }
#'      
#'     Bayesian estimation with \pkg{TreeBUGS}
#'     \itemize{
#'       \item{\code{"simple"}: }{Bayesian estimation, no pooling (C++,
#'         \link[TreeBUGS]{simpleMPT})}
#'       \item{\code{"simple_pooling"}: }{Bayesian estimation, complete pooling 
#'         (C++, \link[TreeBUGS]{simpleMPT})}
#'       \item{\code{"trait"}: }{latent-trait model, partial pooling (JAGS, 
#'         \link[TreeBUGS]{traitMPT})}
#'       \item{\code{"trait_uncorrelated"}: }{latent-trait model without
#'         correlation parameters, partial pooling (JAGS,
#'         \link[TreeBUGS]{traitMPT})}
#'       \item{\code{"beta"}: }{beta-MPT model, partial pooling (JAGS, 
#'         \link[TreeBUGS]{betaMPT})}
#'       \item{\code{"betacpp"}: }{beta-MPT model, partial pooling (C++, 
#'         \link[TreeBUGS]{betaMPTcpp})}
#'     }
#'   }
#'   \subsection{Frequentist/Maximum-Likelihood Methods}{
#'     For the \emph{complete pooling asymptotic approach}, the group-level parameter
#'     estimates and goodness-of-fit statistics are the maximum-likelihood and
#'     G-squared values returned by \code{MPTinR}. The parameter differences are
#'     based on these values, the standard errors of the difference is simply
#'     the pooled standard error of the individual parameters. The overall fit
#'     (column \code{gof}) is based on an additional fit to the completely
#'     aggregated data.
#'     
#'     For the \emph{no pooling asymptotic approach}, the individual-level
#'     maximum-likelihood estimates are reported in column \code{est_indiv} and
#'     \code{gof_indiv} and provide the basis for the other results. The
#'     group-level parameters are simply the means of the individual-level
#'     parameters, the SE is the SE of the mean for these parameter (i.e.,
#'     SD/sqrt(N), where N excludes parameters estimated as NA), and the CI is
#'     based on mean and SE. The group-level and overall fit is the sum of the
#'     individual G-squares, sum of individual-level df, and corresponding
#'     chi-square df. The difference between the conditions and corresponding
#'     statistics are based on a t-test comparing the individual-level
#'     estimates. The CIs of the difference are based on the SEs (which are
#'     derived from an equivalent linear model).
#'     
#'     The individual-level estimates of the \code{bootstrap based no-pooling}
#'     approaches are identical to the asymptotic ones. However, the SE is the
#'     SD of the bootstrapped distribution of parameter estimates, the CIs are
#'     the corresponding quantiles of the bootstrapped distribution, and the
#'     p-value is obtained from the bootstrapped G-square distribution. The
#'     group-level estimates are also the mean of the individual-level
#'     estimates, however, after excluding those parameters that are empirically
#'     not identified based on the CIs derived from the bootstrap distribution.
#'     Specifically, we calculate the range from maximum and minimum CI value
#'     and exclude those individuals for which the range is larger than
#'     \code{mpt_options()$max_ci_indiv} which defaults to \code{0.99}. (PLEASE
#'     CHECK)
#'   }
#' 
#' @return A \code{tibble} with one row per estimation \code{method} and the
#'   following columns:
#' \enumerate{
#'   \item \code{model}: Name of model file (copied from \code{model} argument),
#'   \code{character}
#'   \item \code{dataset}: Name of data set (copied from \code{dataset}
#'   argument), \code{character}
#'   \item \code{pooling}: \code{character} specifying the level of pooling with
#'   three potential values: \code{c("complete", "no", "partial")}
#'   \item \code{package}: \code{character} specifying the package used for
#'   estimation with two potential values: \code{c("MPTinR", "TreeBUGS")}
#'   \item \code{method}: \code{character} specifying the method used with the
#'   following potential values: \code{c("asymptotic", "PB/MLE", "NPB/MLE",
#'   "simple", "trait", "trait_uncorrelated", "beta", "betacpp")}
#'   \item \code{est_group}: Group-level parameter estimates per condition/group.
#'   \item \code{est_indiv}: Individual-level parameter estimates (if provided
#'   by method).
#'   \item \code{test_between}: Parameter differences between the levels of the
#'   between-subjects condition (if specified).
#'   \item \code{gof}: Overall goodness of fit across all individuals.
#'   \item \code{gof_group}: Group-level goodness of fit.
#'   \item \code{gof_indiv}: Individual-level goodness of fit.
#'   \item \code{test_homogeneity}: Chi-square based test of participant
#'   homogeneity proposed by Smith and Batchelder (2008). This test is the same
#'   for each estimation method.
#'   \item \code{convergence}: Convergence information provided by the
#'   respective estimation method. For the asymptotic frequentist methods this
#'   is a \code{tibble} with rank of the Fisher matrix, the number of parameters
#'   (which should match the rank of the Fisgher matrix), and the convergence
#'   code provided by the optimization algorithm (which is
#'   \code{\link{nlminb}}). The boostrap methods contain an additional column,
#'   \code{parameter}, that contains the information which (if any) parameters
#'   are empirically non-identifiable based on the bootstrapped distribution of
#'   parameter estimates (see above for exact description). For the Bayesian
#'   methods this is a \code{tibble} containing information of the posterior
#'   dsitribution (i.e., mean, quantiles, SD, SE, \code{n.eff}, and R-hat) for
#'   each parameter.
#'   \item \code{estimation}: Time it took for each estimation method and group.
#'   \item \code{options}: Options used for estimation. Obtained by running
#'   \code{\link{mpt_options}()}
#' }
#' 
#' With the exception of the first five columns (i.e., after \code{method}) all
#' columns are \code{list} columns typically holding one \code{tibble} per cell.
#' The simplest way to analyze the results is separately per column using
#' \code{\link[tidyr]{unnest}}. Examples for this are given below.
#'
#' @references 
#'   Smith, J. B., & Batchelder, W. H. (2008). Assessing individual differences
#'   in categorical data. \emph{Psychonomic Bulletin & Review}, 15(4), 713-731.
#'   \url{https://doi.org/10.3758/PBR.15.4.713}

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