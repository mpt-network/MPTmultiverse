


mpt_mptinr <- function(
  dataset    # name of data file
  , data     # data.frame
  , model    # name of EQN file
  , method   # analysis approaches to be conducted
  , col_id = "id"
  , col_condition = "condition"
){

  prepared <- prep_data_fitting(data = data,
                            model_file = model,
                            col_id = col_id, 
                            col_condition = col_condition)
  
  res <- list()
  
  if(any(method %in% c("asymptotic_complete"))) {
    res[["complete_pooling"]] <- mpt_mptinr_complete(
      dataset = dataset
      , prepared = prepared
      , model = model
      , method = intersect(method, "asymptotic_complete")
      , col_id = col_id
      , col_condition = col_condition
    )
  }
  
  if(any(method %in% c("asymptotic_no", "pb_no"))) {
    res[["no_pooling"]] <- mpt_mptinr_no(
      dataset = dataset
      , prepared = prepared
      , model = model
      , method = intersect(method, c("asymptotic_no", "pb_no"))
      , col_id = col_id
      , col_condition = col_condition
    )
  }
  
  # return
  dplyr::bind_rows(res)
}


  
################
## no pooling ##
################

#' @importFrom rlang .data sym
#' @importFrom magrittr %>%
#' @keywords internal

mpt_mptinr_no <- function(
  dataset
  , prepared
  , model
  , method
  , col_id
  , col_condition
) {

  OPTIONS <- getOption("MPTmultiverse")
  MPTINR_OPTIONS <- OPTIONS$mptinr
  CI_SIZE <- OPTIONS$ci_size
  MAX_CI_INDIV <- OPTIONS$max_ci_indiv
  
  
  cl <- parallel::makeCluster(rep("localhost", MPTINR_OPTIONS$n.CPU))
  parallel::clusterEvalQ(cl, library("MPTinR"))
  
  no_pooling <- make_results_row(model = model,
                                 dataset = dataset,
                                 pooling = "no",
                                 package = "MPTinR",
                                 method = "PB/MLE",
                                 data = prepared$data,
                                 parameters = prepared$parameters)
  
  no_pooling2 <- make_results_row(model = model,
                                 dataset = dataset,
                                 pooling = "no",
                                 package = "MPTinR",
                                 method = "asymptotic",
                                 data = prepared$data,
                                 parameters = prepared$parameters)
  
  fit_mptinr <- MPTinR::fit.mpt(prepared$data[,prepared$col_freq],
                        model.filename = model,
                        n.optim = MPTINR_OPTIONS$n.optim,
                        fit.aggregated = FALSE,
                        show.messages = FALSE, output = "full", 
                        ci = (1 - stats::pnorm(1))*2*100)
  
  convergence <- data.frame(
    id = prepared$data[, col_id]
    , condition = prepared$data[,col_condition]
    , fit_mptinr$model.info$individual[,1:2]
    , convergence = vapply(fit_mptinr$best.fits$individual, FUN = function(x) x$convergence, 0)
  )
  
  no_pooling$gof_indiv[[1]]$type <- "pb-G2"
  no_pooling$gof_indiv[[1]]$focus <- "mean"
  no_pooling$gof_indiv[[1]]$stat_obs <- fit_mptinr$goodness.of.fit$individual$G.Squared
  no_pooling$gof_indiv[[1]]$stat_df <- fit_mptinr$goodness.of.fit$individual$df
  
  no_pooling2$gof_indiv[[1]]$type <- "G2"
  no_pooling2$gof_indiv[[1]]$focus <- "mean"
  no_pooling2$gof_indiv[[1]]$stat_obs <- fit_mptinr$goodness.of.fit$individual$G.Squared
  no_pooling2$gof_indiv[[1]]$stat_df <- fit_mptinr$goodness.of.fit$individual$df
  no_pooling2$gof_indiv[[1]]$p <- fit_mptinr$goodness.of.fit$individual$p.value

  fit_pb <- parallel::clusterApplyLB(
    cl
    , seq_len(nrow(prepared$data))
    , get_pb_output
    , fit_mptinr = fit_mptinr
    , data = prepared$data
    , model_file = model
    , col_freq = prepared$col_freq
    , MPTINR_OPTIONS = MPTINR_OPTIONS
  )
  
  ## make est_indiv and gof_indiv
  for (i in seq_len(nrow(prepared$data))) {
    
    for (p in prepared$parameters) {
      no_pooling$est_indiv[[1]][
        no_pooling$est_indiv[[1]]$id == prepared$data[i,"id"] &
          no_pooling$est_indiv[[1]]$parameter == p, "est" ] <-
        fit_mptinr$parameters$individual[p,"estimates",i]
      
      no_pooling2$est_indiv[[1]][
        no_pooling2$est_indiv[[1]]$id == prepared$data[i,"id"] &
          no_pooling2$est_indiv[[1]]$parameter == p, "est" ] <-
        fit_mptinr$parameters$individual[p,"estimates",i]
      no_pooling2$est_indiv[[1]][
        no_pooling2$est_indiv[[1]]$id == prepared$data[i,"id"] &
          no_pooling2$est_indiv[[1]]$parameter == p, "se" ] <-
        fit_mptinr$parameters$individual[p, "upper.conf",i] - 
        fit_mptinr$parameters$individual[p,"estimates",i]
      
      no_pooling$est_indiv[[1]][
        no_pooling$est_indiv[[1]]$id == prepared$data[i,"id"] &
          no_pooling$est_indiv[[1]]$parameter == p, prepared$cols_ci ] <-
        stats::quantile(fit_pb[[i]]$parameters$individual[p,"estimates",], probs = CI_SIZE)
      no_pooling$est_indiv[[1]][
        no_pooling$est_indiv[[1]]$id == prepared$data[i,"id"] &
          no_pooling$est_indiv[[1]]$parameter == p, "se" ] <-
        stats::sd(fit_pb[[i]]$parameters$individual[p,"estimates",]) 
    }
    
    # gof_indiv
    no_pooling$gof_indiv[[1]][
      no_pooling$gof_indiv[[1]]$id == prepared$data[i,"id"], "p" ] <-
      (sum(fit_pb[[i]]$goodness.of.fit$individual$G.Squared >
             fit_mptinr$goodness.of.fit$individual[i,"G.Squared"]) + 1) /
      (MPTINR_OPTIONS$bootstrap_samples + 1)
    
  }
  
  for (i in seq_along(CI_SIZE)) {
    no_pooling2$est_indiv[[1]][, prepared$cols_ci[i]] <-
      no_pooling2$est_indiv[[1]][,"est"] +
      stats::qnorm(CI_SIZE[i])*no_pooling2$est_indiv[[1]][,"se"]
  }
  
  #### make est_group ####
  
  tmp <- no_pooling$est_indiv[[1]]
  tmp$range_ci <- tmp[, prepared$cols_ci[length(prepared$cols_ci)]][[1]] - 
    tmp[, prepared$cols_ci[1]][[1]]
  
  non_identified_pars <-  tmp %>%
    dplyr::filter(.data$range_ci > MAX_CI_INDIV) %>% 
    dplyr::group_by(.data$id) %>% 
    dplyr::summarise(parameter = paste0(.data$parameter, collapse = ", ")) %>% 
    dplyr::ungroup()
  
  no_pooling$convergence <- 
    list(tibble::as_tibble(dplyr::left_join(convergence, non_identified_pars, by = "id")))
  
  if (nrow(non_identified_pars) > 0) {
    warning("MPTinR-no: IDs and parameters with PB-CIs > ",
            MAX_CI_INDIV, " (i.e., non-identified):\n", 
            apply(non_identified_pars, 
                  1, function(x) paste0(x[["id"]], ": ", x[["parameter"]], "\n") ),
            call. = FALSE)    
  }

  est_group <- tmp %>%
    dplyr::filter(.data$range_ci < MAX_CI_INDIV) %>%
    dplyr::group_by(.data$condition, .data$parameter) %>%
    dplyr::summarise(
      estN = mean(.data$est)
      , se = stats::sd(.data$est) / sqrt(sum(!is.na(.data$est)))
      , quant = list(as.data.frame(t(stats::quantile(.data$est, prob = CI_SIZE))))
    ) %>%
    tidyr::unnest(.data$quant) %>%
    dplyr::ungroup() %>%
    dplyr::rename(est = .data$estN)
  colnames(est_group)[
    (length(colnames(est_group))-length(CI_SIZE)+1):length(colnames(est_group))
    ] <- prepared$cols_ci
  
  no_pooling$est_group[[1]] <-
    dplyr::right_join(est_group,
               no_pooling$est_group[[1]][, c("condition", "parameter")],
               by = c("condition", "parameter"))
  
  
    est_group2 <- no_pooling2$est_indiv[[1]] %>%
    dplyr::group_by(.data$condition, .data$parameter) %>%
    dplyr::summarise(estN = mean(.data$est),
              se = stats::sd(.data$est) / sqrt(sum(!is.na(.data$est))),
              quant = list(as.data.frame(t(stats::quantile(.data$est, prob = CI_SIZE))))) %>%
    tidyr::unnest(.data$quant) %>%
    dplyr::ungroup() %>%
    dplyr::rename(est = .data$estN)
  colnames(est_group2)[
    (length(colnames(est_group2))-length(CI_SIZE)+1):length(colnames(est_group2))
    ] <- prepared$cols_ci
  
  no_pooling2$est_group[[1]] <-
    dplyr::right_join(est_group2,
               no_pooling2$est_group[[1]][,c("condition", "parameter")],
               by = c("condition", "parameter"))
  
  #### make gof_group ####
  no_pooling$gof_group[[1]]$type <- "pb-G2"
  no_pooling$gof_group[[1]]$focus <- "mean"
  
  tmp <- fit_mptinr$goodness.of.fit$individual
  tmp$condition <- as.factor(prepared$data$condition)
  gof_group <- tmp %>%
    dplyr::group_by(.data$condition) %>%
    dplyr::summarise(stat_obs = sum(.data$G.Squared),
              stat_df = sum(.data$df))
  gof_group$p <- NA_real_
  
  g2_all <- vapply(fit_pb,
                   function(x) x$goodness.of.fit$individual$G.Squared,
                   rep(0, MPTINR_OPTIONS$bootstrap_samples))
  
  g2_cond <- vector("list", length(prepared$conditions))
  
  for (i in seq_along(prepared$conditions)) {
    g2_cond[[i]] <- apply(g2_all[ , 
                                 prepared$data$condition == 
                                   prepared$conditions[i]], 1, sum)
    no_pooling$gof_group[[1]][ 
      no_pooling$gof_group[[1]]$condition == 
        prepared$conditions[i], "stat_obs" ] <- 
      gof_group[ gof_group$condition == prepared$conditions[i], "stat_obs"]
    no_pooling$gof_group[[1]][ no_pooling$gof_group[[1]]$condition == 
                                 prepared$conditions[i], "stat_df" ] <- 
      gof_group[ gof_group$condition == prepared$conditions[i], "stat_df"]
    no_pooling$gof_group[[1]][ no_pooling$gof_group[[1]]$condition == 
                                 prepared$conditions[i], "p" ] <-
      (sum(gof_group[ gof_group$condition == 
                        prepared$conditions[i], "stat_obs"][[1]] <
             g2_cond[[i]]) + 1) / (MPTINR_OPTIONS$bootstrap_samples + 1)
  }
  
  #### make gof_group2 ####
  no_pooling2$gof_group[[1]]$type <- "G2"
  no_pooling2$gof_group[[1]]$focus <- "mean"
  
  gof_group2 <- tmp %>%
    dplyr::group_by(.data$condition) %>%
    dplyr::summarise(stat_obs = sum(.data$G.Squared),
              stat_pred = NA_real_,
              stat_df = sum(.data$df))
  gof_group2$p <- stats::pchisq(q = gof_group2$stat_obs, 
                         df = gof_group2$stat_df, 
                         lower.tail = FALSE)
  gof_group2$condition <- factor(gof_group2$condition)
  no_pooling2$gof_group[[1]] <- 
    dplyr::right_join(no_pooling2$gof_group[[1]][, c("condition", "type", "focus")],
             gof_group2,
               by = c("condition"))
  
  no_pooling2$gof[[1]]$type <- "G2"
  no_pooling2$gof[[1]]$focus <- "mean"
  no_pooling2$gof[[1]]$stat_obs <- fit_mptinr$goodness.of.fit$sum$G.Squared
  no_pooling2$gof[[1]]$stat_df <- fit_mptinr$goodness.of.fit$sum$df
  no_pooling2$gof[[1]]$p <- stats::pchisq(q = no_pooling2$gof[[1]]$stat_obs, 
                                   df = no_pooling2$gof[[1]]$stat_df, 
                                   lower.tail = FALSE)
  
  #### make gof ####
  no_pooling$gof[[1]]$type <- "pb-G2"
  no_pooling$gof[[1]]$focus <- "mean"
  no_pooling$gof[[1]]$stat_obs <- fit_mptinr$goodness.of.fit$sum$G.Squared
  no_pooling$gof[[1]]$stat_df <- fit_mptinr$goodness.of.fit$sum$df
  
  g2_all <- vapply(fit_pb,
                   function(x) x$goodness.of.fit$individual$G.Squared,
                   rep(0, MPTINR_OPTIONS$bootstrap_samples))
  
  g2_cond <- apply(g2_all, 1, sum)
  no_pooling$gof[[1]]$p <-
    (sum(no_pooling$gof[[1]]$stat_obs < g2_cond) + 1) /
    (MPTINR_OPTIONS$bootstrap_samples + 1)
  
  ### test between ###
  
  for (i in seq_len(nrow(no_pooling$test_between[[1]]))) {
    tmp_par <- no_pooling$test_between[[1]]$parameter[i]
    tmp_c1 <- no_pooling$test_between[[1]]$condition1[i]
    tmp_c2 <- no_pooling$test_between[[1]]$condition2[i]
    
    tmp_df <- droplevels(no_pooling$est_indiv[[1]][ 
      no_pooling$est_indiv[[1]]$parameter == tmp_par & 
        no_pooling$est_indiv[[1]]$condition %in% 
        c(as.character(tmp_c1), as.character(tmp_c2)) , ])
    
    tmp_t <- stats::t.test(tmp_df[ tmp_df$condition == tmp_c1,  ]$est, 
                    tmp_df[ tmp_df$condition == tmp_c2,  ]$est)
    
    tmp_lm <- stats::lm(est ~ condition, tmp_df)
    
    tmp_se <- stats::coef(stats::summary.lm(tmp_lm))[2,"Std. Error"]
    
    no_pooling$test_between[[1]][ 
      no_pooling$test_between[[1]]$parameter == tmp_par , 
      c("est_diff" , "se", "p") ] <- 
      c(diff(rev(tmp_t$estimate)), tmp_se, tmp_t$p.value)
    
    no_pooling$test_between[[1]][ 
      no_pooling$test_between[[1]]$parameter == tmp_par, prepared$cols_ci] <- 
      no_pooling$test_between[[1]][ 
        no_pooling$test_between[[1]]$parameter == tmp_par ,]$est_diff + 
      stats::qnorm(CI_SIZE)* no_pooling$test_between[[1]][ 
        no_pooling$test_between[[1]]$parameter == tmp_par ,]$se
    
  }
  
  ### copy information
  no_pooling2$convergence <- no_pooling$convergence
  no_pooling2$test_between <- no_pooling$test_between
  
  parallel::stopCluster(cl)
  
  return(dplyr::bind_rows(no_pooling, no_pooling2))
}


#' Parametric Bootstrap for MPT
#'
#' Helper function for creating parametric-bootstrap confidence intervals.
#'
#' @keywords internal

get_pb_output <- function(
  i
  , fit_mptinr
  , data
  , model_file
  , col_freq
  , MPTINR_OPTIONS) {
  
  gen_data <- MPTinR::gen.data(
    fit_mptinr$parameters$individual[,"estimates", i]
    , samples = MPTINR_OPTIONS$bootstrap_samples
    , model.filename = model_file
    , data = unlist(data[i, col_freq])
  )
  MPTinR::fit.mpt(gen_data,
            model.filename = model_file,
            fit.aggregated = FALSE,
            n.optim = MPTINR_OPTIONS$n.optim,
            show.messages = FALSE)
}


######################
## complete pooling ##
######################

#' @importFrom magrittr %>%
#' @keywords internal

mpt_mptinr_complete <- function(dataset, 
                                prepared, 
                                model,
                                method,
                                col_id, 
                                col_condition) {
  OPTIONS <- getOption("MPTmultiverse")
  MPTINR_OPTIONS <- OPTIONS$mptinr
  CI_SIZE <- OPTIONS$ci_size
  MAX_CI_INDIV <- OPTIONS$max_ci_indiv
  
  
  complete_pooling <- make_results_row(model = model,
                                       dataset = dataset,
                                       pooling = "complete",
                                       package = "MPTinR",
                                       method = "asymptotic",
                                       data = prepared$data,
                                       parameters = prepared$parameters)
  
  complete_pooling$est_indiv <- list(tibble::tibble())
  complete_pooling$gof_indiv <- list(tibble::tibble())
  
  #### fully aggregated:
  
  fit_mptinr_agg <- MPTinR::fit.mpt(colSums(prepared$data[,prepared$col_freq]),
                            model.filename = model,
                            n.optim = MPTINR_OPTIONS$n.optim,
                            show.messages = FALSE, output = "full")
  
  ## gof
  
  complete_pooling$gof[[1]][1,"type"] <- "G2"
  complete_pooling$gof[[1]][1,"focus"] <- "mean"
  complete_pooling$gof[[1]][1,"stat_obs"] <-
    fit_mptinr_agg$goodness.of.fit$G.Squared
  complete_pooling$gof[[1]][1,"stat_df"] <-
    fit_mptinr_agg$goodness.of.fit$df
  complete_pooling$gof[[1]][1,"p"] <-
    fit_mptinr_agg$goodness.of.fit$p
  
  #### aggregated by condition
  
  complete_pooling$gof_group[[1]][,"type"] <- "G2"
  complete_pooling$gof_group[[1]][,"focus"] <- "mean"
  
  convergence <- vector("list", 1 + length(prepared$conditions))
  names(convergence) <- c("aggregated", prepared$conditions)
  convergence$aggregated <- tibble::as_tibble(data.frame(
             fit_mptinr_agg$model.info[,1:2],
             convergence = fit_mptinr_agg$best.fits[[1]]$convergence))
  
  
  for (i in seq_along(prepared$conditions)) {
    fit_mptinr_tmp <- MPTinR::fit.mpt(colSums(
      prepared$freq_list[[prepared$conditions[i]]][,prepared$col_freq]),
                              model.filename = model,
                              n.optim = MPTINR_OPTIONS$n.optim,
                              show.messages = FALSE,
                              output = "full")
    
    
    complete_pooling$gof_group[[1]][
      complete_pooling$gof_group[[1]]$condition == 
        prepared$conditions[i] ,"stat_obs"] <-
      fit_mptinr_tmp$goodness.of.fit$G.Squared
    complete_pooling$gof_group[[1]][
      complete_pooling$gof_group[[1]]$condition == 
        prepared$conditions[i], "stat_df"] <-
      fit_mptinr_tmp$goodness.of.fit$df
    complete_pooling$gof_group[[1]][
      complete_pooling$gof_group[[1]]$condition == 
        prepared$conditions[i], "p"] <-
      fit_mptinr_tmp$goodness.of.fit$p
    
    complete_pooling$est_group[[1]][
      complete_pooling$est_group[[1]]$condition == 
        prepared$conditions[i], "est"
      ] <- fit_mptinr_tmp$parameters[
        complete_pooling$est_group[[1]][
          complete_pooling$est_group[[1]]$condition == 
            prepared$conditions[i],
          ]$parameter, "estimates" ]
    
    par_se <- sqrt(diag(solve(fit_mptinr_tmp$hessian[[1]])))
    names(par_se) <- rownames(fit_mptinr_tmp$parameters)
    
    complete_pooling$est_group[[1]][
      complete_pooling$est_group[[1]]$condition == 
        prepared$conditions[i], "se"
      ] <- par_se[
        complete_pooling$est_group[[1]][
          complete_pooling$est_group[[1]]$condition == 
            prepared$conditions[i],
          ]$parameter]
    
    convergence[[prepared$conditions[i]]] <- tibble::as_tibble(data.frame(
             fit_mptinr_tmp$model.info[,1:2],
             convergence = fit_mptinr_tmp$best.fits[[1]]$convergence))
  }
  
  for (i in seq_along(CI_SIZE)) {
    complete_pooling$est_group[[1]][, prepared$cols_ci[i]] <-
      complete_pooling$est_group[[1]][,"est"] +
      stats::qnorm(CI_SIZE[i])*complete_pooling$est_group[[1]][,"se"]
  }
  
  ### test between ###
  tmp_pars <- complete_pooling$est_group[[1]]
  
  for (i in seq_len(nrow(complete_pooling$test_between[[1]]))) {
    
    tmp_par <- complete_pooling$test_between[[1]]$parameter[i]
    tmp_c1 <- complete_pooling$test_between[[1]]$condition1[i]
    tmp_c2 <- complete_pooling$test_between[[1]]$condition2[i]
  
    complete_pooling$test_between[[1]][i, "est_diff"] <- 
      tmp_pars[tmp_pars$condition == tmp_c1 & 
               tmp_pars$parameter == tmp_par, ]$est -
      tmp_pars[tmp_pars$condition == tmp_c2 & 
               tmp_pars$parameter == tmp_par, ]$est
    
  }
  
  tmp <- names(convergence)
  convergence <- do.call("rbind", convergence)
  convergence <- dplyr::bind_cols(condition = factor(tmp), 
                           convergence)
  
  complete_pooling$convergence <- list(convergence)
  warn_conv <- convergence$convergence != 0
  if (any(warn_conv)) {
    warning("MPTinR-complete: Convergence code != 0 for: ", 
            paste0(names(warn_conv)[warn_conv], collapse = ", "), 
            call. = FALSE)
  }
  
  return(complete_pooling)
}
