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
  , id = NULL
  , condition = NULL
) {
  
  OPTIONS <- getOption("MPTmultiverse")
  
  CI_SIZE <- OPTIONS$ci_size
  MAX_CI_INDIV <- OPTIONS$max_ci_indiv
  
  prepared <- MPTmultiverse:::prep_data_fitting(
    data = data
    , model_file = model
    , id = id
    , condition = condition
  )
  
  results_row <- make_results_row(
    model = model
    , dataset = dataset
    , pooling = "partial"
    , package = "HMMTreeR"
    , method = "latent_class"
    , data = prepared$data
    , parameters = prepared$parameters
    , id = id
    , condition = condition
  )
  
  # Currently, no individual parameter estimates (and model fits) are calculated.
  # Therefore, these should be empty tibbles.
  results_row$est_indiv <- list(tibble::tibble())
  results_row$gof_indiv <- list(tibble::tibble())
  
  
  # ----------------------------------------------------------------------------
  # Aggregate analyses: Ignore between-subjects condition
  
  # create a temporary directory
  tmp_dir_name <- paste(c("HMMTreeR-tmp-directory-", sample(c(letters, LETTERS), 20, replace = TRUE)), collapse = "")
  dir.create(tmp_dir_name)
  
  # copy .eqn file to tmp dir
  file.copy(from = model, to = tmp_dir_name, copy.mode = FALSE)
  simplify_eqn(model_filename = model, eqn_filename = file.path(tmp_dir_name, model), data = data, id = id, condition = condition)
  
  # write data of condition group to tab-separated file
  data_file <- file.path(tmp_dir_name, dataset)
  id_col <- data.frame(id = rep(dataset, nrow(data)), stringsAsFactors = FALSE)
  write.table(file = data_file, x = cbind(id_col, data[, setdiff(colnames(data), c(id, condition))]), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  print(data_file)
  print(file.path(tmp_dir_name, model))
  
  
  t0 <- Sys.time()
  res <- HMMTreeR::lc(
    model = file.path(tmp_dir_name, model)
    , data = data_file
    , nsubj = nrow(data)
    , max_classes = getOption("MPTmultiverse")$hmmtree$max_classes
    , runs = getOption("MPTmultiverse")$hmmtree$n.optim
    , fisher_information = getOption("MPTmultiverse")$hmmtree$fisher_information
  )
  t1 <- as.numeric(Sys.time() - t0)
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
  estimation_time <- list()
  
  for (j in prepared$conditions) {
    # create a temporary directory
    tmp_dir_name <- paste(c("HMMTreeR-tmp-directory-", sample(c(letters, LETTERS), 20, replace = TRUE)), collapse = "")
    dir.create(tmp_dir_name)
    
    # copy .eqn file to tmp dir
    simplify_eqn(model_filename = model, eqn_filename = file.path(tmp_dir_name, model), data = data, id = id, condition = condition)
    
    # write data of condition group to tab-separated file
    data_file <- file.path(tmp_dir_name, dataset)
    id_col <- data.frame(id = rep(dataset, nrow(prepared$freq_list[[j]])))
    write.table(file = data_file, x = cbind(id_col, prepared$freq_list[[j]]), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    t0 <- Sys.time()
    res <- HMMTreeR::lc(
      model = file.path(tmp_dir_name, model)
      , data = data_file
      , nsubj = nrow(prepared$freq_list[[j]])
      , max_classes = getOption("MPTmultiverse")$hmmtree$max_classes
      , runs = getOption("MPTmultiverse")$hmmtree$n.optim
      , fisher_information = getOption("MPTmultiverse")$hmmtree$fisher_information
    )
    estimation_time[[j]] <- as.numeric(Sys.time() - t0)
    
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
  
  est_group <- dplyr::bind_rows(est_group)
  gof_group <- dplyr::bind_rows(gof_group)
  
  results_row$est_group[[1]] <- est_group
  results_row$gof_group[[1]] <- gof_group
  

  # test_between ----
  test_between <- results_row$test_between[[1]]
  
  for (i in seq_len(nrow(test_between))) {
    
    p <- test_between$parameter[i]
    c1 <- test_between$condition1[i]
    c2 <- test_between$condition2[i]
    
    test_between[i, "est_diff"] <- 
      est_group[est_group$condition == c1 & est_group$parameter == p, ]$est -
      est_group[est_group$condition == c2 & est_group$parameter == p, ]$est
    
    test_between[i, "se"] <- sqrt(
      (est_group[est_group$condition == c1 & est_group$parameter == p, ]$se)^2 +
      (est_group[est_group$condition == c2 & est_group$parameter == p, ]$se)^2
    )
  }
  
  for (k in OPTIONS$ci_size) {
    test_between[[paste0("ci_", k)]] <- test_between$est_diff + test_between$se * qnorm(p = k)
  }

  results_row$test_between[[1]] <- test_between

  estimation_time <- unlist(estimation_time)
  
  results_row$estimation[[1]] <- tibble::tibble(
    condition = c("complete_data", names(estimation_time))
    , time_difference = c(t1, unname(estimation_time))
  )
  
  # return ----
  results_row
}

#' @keywords internal

simplify_eqn <- function(model_filename, eqn_filename, data = data, id, condition) {
  read_lines <- readLines(model_filename, warn = FALSE)
  n_terms <- as.integer(read_lines[[1]])
  if(n_terms!=(length(read_lines)-1)) stop("eqn file seems to be corrupted")
  repeat{
    read_lines <- gsub(read_lines, pattern = "  |\t", replacement = " ")
    if(!any(grepl(read_lines, pattern = "  |\t"))) {
      break
    }
  }
  read_lines <- sub(read_lines, pattern = " ", replacement = ";")
  read_lines <- sub(read_lines, pattern = " ", replacement = ";")
  read_lines <- strsplit(read_lines[-1], split = ";")
  
  model <- data.frame(
    tree = unlist(lapply(X = read_lines, FUN = function(x){x[1]}))
    , cat = unlist(lapply(X = read_lines, FUN = function(x){x[2]}))
    , term = unlist(lapply(X = read_lines, FUN = function(x){x[3]}))
    , stringsAsFactors = FALSE
  )
  cat_cols <- setdiff(colnames(data), c(id, condition))
  cat_id <- seq_along(cat_cols)
  names(cat_id) <- cat_cols
  
  for(i in 1:length(unique(model$tree))) {
    model$tree[model$tree==unique(model$tree)[i]] <- i
  }

  model$cat <- cat_id[model$cat]
  writeLines(
    text = paste(c(n_terms, paste(model$tree, model$cat, model$term, sep = " ")), collapse = "\n")
    , con = eqn_filename
    , sep=""
  )

}
