#' Fit latent-class multinomial models
#'
#' Does a lot of nice stuff
#'
#' @param dataset Character.
#' @param data A \code{data.frame} containing the data.
#' @param model A model definition, typically the path to an \code{.eqn} file.
#' @param id Character. Name of the column that contains the subject identifier.
#'   If not specified, it is assumed that each row represents observations from one participant.
#' @param condition Character. Name of the column specifying a between-subjects factor.
#'   If not specified, no between-subjects comparisons are performed.
#'
#' @keywords internal


fit_lc <- function(
  dataset
  , data
  , model
  , col_id = NULL
  , col_condition = NULL
) {
  
  OPTIONS <- getOption("MPTmultiverse")
  
  CI_SIZE <- OPTIONS$ci_size
  MAX_CI_INDIV <- OPTIONS$max_ci_indiv
  
  prepared <- MPTmultiverse:::prep_data_fitting(
    data = data
    , model_file = model
    , col_id = col_id
    , col_condition = col_condition
  )
  
  results_row <- MPTmultiverse:::make_results_row(
    model = model
    , dataset = dataset
    , pooling = "partial"
    , package = "HMMTreeR"
    , method = "latent_class"
    , data = prepared$data
    , parameters = prepared$parameters
  )
  
  # Currently, no individual parameter estimates (and model fits) are calculated.
  # Therefore, these should be empty tibbles.
  results_row$est_indiv <- list(tibble::tibble())
  results_row$gof_indiv <- list(tibble::tibble())
  
  
  # ----------------------------------------------------------------------------
  # Aggregate analyses: Ignore between-subjects condition
    
  res <- HMMTreeR::lc(
    model = model
    , data = data_file
    , nsubj = nrow(data)
    , max_classes = getOption("MPTmultiverse")$hmmtree$max_classes
    , nruns = getOption("MPTmultiverse")$hmmtree$n.optim
    , fisher_information = getOption("MPTmultiverse")$hmmtree$fisher_information
  )
  
  fit_stats <- HMMTreeR::fit_statistics(res[[length(res)]]) # choose winning model
  
  required_stats <- c("M1", "M2", "S1", "S2")
  
  gof <- tibble::tibble(
    type = required_stats
    , focus = rep(c("mean", "cov"), each = 2)
    , stat_obs = as.numeric(fit_stats[, required_stats])
    , stat_pred = NA_real_
    , stat_df = as.numeric(fit_stats[, paste0("df_", required_stats)])
    , p = ifelse(stat_df <= 0, NA_real_, pchisq(q = stat_obs, df = stat_df))
  )
  
  results_row$gof[[1]] <- gof
  
  # ----------------------------------------------------------------------------
  # Analyses by between-subjects condition
  
  est_group <- list()
  gof_group <- list()
  
  for (j in prepared$conditions) {
    # create a temporary directory
    tmp_dir_name <- paste(c("HMMTreeR-tmp-directory-", sample(c(letters, LETTERS), 20, replace = TRUE)), collapse = "")
    dir.create(tmp_dir_name)
    
    # copy .eqn file to tmp dir
    file.copy(from = model, to = tmp_dir_name, copy.mode = FALSE)
    
    # write data of condition group to tab-separated file
    data_file <- file.path(tmp_dir_name, dataset)
    write.table(file = data_file, x = prepared$freq_list[[j]], sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    res <- HMMTreeR::lc(
      model = model
      , data = data_file
      , nsubj = nrow(data)
      , max_classes = getOption("MPTmultiverse")$hmmtree$max_classes
      , nruns = getOption("MPTmultiverse")$hmmtree$n.optim
      , fisher_information = getOption("MPTmultiverse")$hmmtree$fisher_information
    )
    
    # parameter estimates ----
    estimates <- HMMTreeR::weighted_means(res[[length(res)]])
    
    est_group[[j]] <- tibble::tibble(
      condition = j
      , parameter = estimates$parameter
      , est = estimates$estimate
      , se = (estimates$upper - estimates$estimate) / qnorm(p = .975)
    )
    
    for (k in OPTIONS$ci_size) {
      est_group[[j]][[paste0("ci_", k)]] <- est_group[[j]]$est + est_group[[j]]$se * qnorm(p = k)
    }
    
    # Overwrite with exact values from HMMTree output, if possible:
    if(.025 %in% OPTIONS$ci_size) {
      est_group[[j]][["ci_0.025"]] <- estimates$lower
    }
    if(.975 %in% OPTIONS$ci_size) {
      est_group[[j]][["ci_0.975"]] <- estimates$upper
    }
    
    
    
    # goodness-of-fit ----
    fit_stats <- HMMTreeR::fit_statistics(res[[length(res)]]) # choose winning model
    
    required_stats <- c("M1", "M2", "S1", "S2")
    
    gof_group[[j]] <- tibble::tibble(
      condition = j
      , type = required_stats
      , focus = rep(c("mean", "cov"), each = 2)
      , stat_obs = as.numeric(fit_stats[, required_stats])
      , stat_pred = NA_real_
      , stat_df = as.numeric(fit_stats[, paste0("df_", required_stats)])
      , p = ifelse(stat_df <= 0, NA_real_, pchisq(q = stat_obs, df = stat_df))
    )
    
  }
  
  results_row$est_group[[1]] <- dplyr::bind_rows(est_group)
  results_row$gof_group[[1]] <- dplyr::bind_rows(gof_group)
  
  # return
  results_row
}
