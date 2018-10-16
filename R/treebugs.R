
# overall goodness of fit across between-subject conditions
aggregate_ppp <- function(ppp_list, stat = "T1"){
  obs <- vapply(ppp_list, "[[", paste0(stat, ".obs"),
                FUN.VALUE = ppp_list[[1]]$T1.obs)
  pred <- vapply(ppp_list, "[[", paste0(stat, ".pred"),
                 FUN.VALUE = ppp_list[[1]]$T1.obs)
  s.obs <- rowSums(obs)
  s.pred <- rowSums(pred)
  c(stat_obs = mean(s.obs), stat_pred = mean(s.pred),
    stat_df = NA, p = mean(s.obs < s.pred))
}

#' @importFrom rlang .data
#' @importFrom magrittr %>%
 
mpt_treebugs <- function (
  method
  , dataset
  , data
  , model
  , id = "id"
  , condition = "condition"
  , core = NULL
){
  all_options <- getOption("MPTmultiverse")
  
  TREEBUGS_MCMC <- all_options$treebugs
  CI_SIZE <- all_options$ci_size
  
  # dlist <- prepare_data(model, data, id = "id", condition = "condition")
  conditions <- levels(factor(data[[condition]]))
  parameters <- MPTinR::check.mpt(model)$parameters
  col_freq <- get_eqn_categories(model)

  data$id <- data[, id]
  data$condition <- data[, condition]
  
  freq_list <- split(data[, col_freq], f = data[, condition])
  pooling <- switch(method, 
                    "simple" = "no", 
                    "simple_pooling" = "complete",
                    "partial")
  
  result_row <- make_results_row(model = model,
                                 dataset = dataset,
                                 pooling = pooling,
                                 package = "TreeBUGS",
                                 method = sub("_pooling","", method, fixed = TRUE),
                                 data = data,
                                 parameters = parameters,
                                 id = id,
                                 condition = condition,
                                 core = core)
  if (method == "simple_pooling"){
    method <- "simple"
    
    # pooling: aggregate across participants
    data <- stats::aggregate(data[, col_freq], list(condition = data$condition), sum)
    data[[condition]] <- data$condition
    data[[id]] <- data$id <- 1:nrow(data)
    if(condition!="condition"){
      data$condition <- NULL
    }
    
    freq_list <- lapply(freq_list, function(x) as.matrix(colSums(x)))
  } 
  if (method == "trait_uncorrelated"){
    method <- "trait"
    prior_args <- list(df = 1, V = NA, xi = "dnorm(0,1)")
  } else {
    prior_args <- NULL
  }
  
  gof_group <- list()
  treebugs_fit <- list()
  estimation_time <- list()
  
  
  for (i in seq_along(conditions)){
    cond <- factor(conditions[i], conditions)
    sel_condition <- data[[condition]] == conditions[i]
    data_group <- data[sel_condition, col_freq]   #freq_list[[i]]
    rownames(data_group) <- data[[id]][sel_condition]
    
    fit_args <- list(eqnfile=model,
                     data = data_group,
                     n.chains = TREEBUGS_MCMC$n.chains,
                     n.iter = TREEBUGS_MCMC$n.iter,
                     n.adapt = TREEBUGS_MCMC$n.adapt,
                     n.burnin = TREEBUGS_MCMC$n.burnin,
                     n.thin = TREEBUGS_MCMC$n.thin)
    if (method %in% c("simple", "betacpp")){
      fit_args["n.adapt"] <- NULL
      fit_args <- c(fit_args, cores = unname(all_options$n.CPU))
    }
    # print(c(fit_args, prior_args))
    t0 <- Sys.time()
    treebugs_function <- ifelse(method == "betacpp", 
                                "TreeBUGS::betaMPTcpp",
                                paste0("TreeBUGS::", method, "MPT"))
    treebugs_fit[[i]] <- do.call(eval(parse(text = treebugs_function)), 
                                 args = c(fit_args, prior_args))
    summ <- treebugs_fit[[i]]$mcmc.summ
    
    # continue MCMC sampling (only for betaMPT and traitMPT)
    ext_cnt <- 0
    try({
      while (
        ext_cnt < TREEBUGS_MCMC$extend_max && method %in% c("beta", "trait") &&
        (any(stats::na.omit(summ[,"Rhat"]) > TREEBUGS_MCMC$Rhat_max)  ||
         any(summ[summ[,"n.eff"] > 0,"n.eff"] < TREEBUGS_MCMC$Neff_min, na.rm = TRUE)) ){
        cat("Drawing additional samples for method = ", method, 
            ". max(Rhat) = ", round(max(stats::na.omit(summ[summ[,"Rhat"] > 0,"Rhat"])), 2),
            " ; min(n.eff) = ", round(min(summ[summ[,"n.eff"] > 0,"n.eff"], na.rm = TRUE), 1), "\n")
        
        treebugs_fit[[i]] <- TreeBUGS::extendMPT(treebugs_fit[[i]],
                                       n.iter = TREEBUGS_MCMC$n.iter,
                                       n.thin = TREEBUGS_MCMC$n.thin,
                                       n.adapt = TREEBUGS_MCMC$n.adapt)
        summ <- treebugs_fit[[i]]$mcmc.summ
        ext_cnt <- ext_cnt + 1
      }
    })
    
    estimation_time[[conditions[i]]] <- Sys.time() - t0
    
    # convergence summary (n.eff / Rhat / all estimates)
    tsum <- tibble::as_tibble(summ) %>% 
      dplyr::mutate(parameter = rownames(summ),
             condition = as.character(cond)) %>% 
      dplyr::select(.data$condition, .data$parameter, .data$Mean : .data$Rhat)
    result_row$convergence[[1]] <- dplyr::bind_rows(result_row$convergence[[1]], tsum)
    
    # parameter estimates
    summMPT <- TreeBUGS::summarizeMPT(treebugs_fit[[i]]$runjags$mcmc,
                                      mptInfo = treebugs_fit[[i]]$mptInfo,
                                      probs = CI_SIZE,
                                      summ = treebugs_fit[[i]]$mcmc.summ)
    
    sel_group <- result_row$est_group[[1]]$condition == conditions[i]
    result_row$est_group[[1]][sel_group,-(1:3)] <-
      summMPT$groupParameters$mean[paste0("mean_", parameters),1:6]
    
    if (pooling != "complete"){
      # # old: array filled into data frame
      # result_row$est_indiv[[1]][sel_ind,-(1:4)] <-
      #   summMPT$individParameters[parameters,,1:(2+length(CI_SIZE))]
      sel_ind <- result_row$est_indiv[[1]]$condition == conditions[i]
      dimnames(summMPT$individParameters)$ID <- rownames(data_group)
      tmp <- summMPT$individParameters[parameters,,1:(2+length(CI_SIZE)), drop = FALSE] %>%
        reshape2::melt() %>% 
        tidyr::spread("Statistic", "value")
      colnames(tmp) <- c("parameter", "id", colnames(result_row$est_indiv[[1]])[-(1:4)])
      tmp[[condition]] <- cond
      result_row$est_indiv[[1]][sel_ind,] <-
        dplyr::left_join(result_row$est_indiv[[1]][sel_ind,] %>%
                    dplyr::select("id", "condition", "parameter", "core"),
                  tmp, by = c("parameter", "id", condition = condition))
    }
    
    gof_group[[i]] <- TreeBUGS::PPP(treebugs_fit[[i]], M = TREEBUGS_MCMC$n.PPP, type = "G2",
                          T2 = pooling != "complete", nCPU = all_options$n.CPU)

    sel_gof <- result_row$gof_group[[1]]$condition == conditions[i]
    result_row$gof_group[[1]][sel_gof, ] <-
      result_row$gof_group[[1]] %>%
      dplyr::filter(.data$condition == conditions[i]) %>%
      dplyr::mutate(
        condition = conditions[i],
        type = "T1",
        focus = "mean",
        stat_obs = mean(gof_group[[i]]$T1.obs),
        stat_pred = mean(gof_group[[i]]$T1.pred),
        p = gof_group[[i]]$T1.p
      )
    
    if (pooling != "complete"){
      result_row$gof_group[[1]] <- tibble::add_row(result_row$gof_group[[1]],
                                           condition = cond,
                                           type = "T2", focus = "cov",
                                           stat_obs = mean(gof_group[[i]]$T2.obs),
                                           stat_pred = mean(gof_group[[i]]$T2.pred),
                                           p = gof_group[[i]]$T2.p)
      
      sel_fog_ind <- result_row$gof_indiv[[1]]$condition == conditions[i]
      result_row$gof_indiv[[1]][sel_fog_ind,] <-
        result_row$gof_indiv[[1]][sel_fog_ind,] %>%
        dplyr::mutate(
          condition = conditions[i],
          type = "T1_G2",
          focus = "mean",
          stat_obs = colMeans(gof_group[[i]]$ind.T1.obs),
          stat_pred = colMeans(gof_group[[i]]$ind.T1.pred),
          p = gof_group[[i]]$ind.T1.p)
    }
  }
  
  # between  subject comparisons
  if (length(conditions) > 1){
    for (i in 1:(length(conditions) - 1)){
      for (j in 2:length(conditions)){
        for(p in parameters){
          test_between <- TreeBUGS::betweenSubjectMPT(treebugs_fit[[i]], treebugs_fit[[j]], 
                                            par1 = p, stat = "x-y")
          test_summ <- TreeBUGS::summarizeMCMC(test_between$mcmc, probs = CI_SIZE)
          bayesp <- mean(do.call("rbind", test_between$mcmc) <= 0)
          
          sel_row <- 
            result_row$test_between[[1]]$parameter == p &
            result_row$test_between[[1]]$condition1 == conditions[i] &
            result_row$test_between[[1]]$condition2 == conditions[j]
          
          result_row$test_between[[1]][sel_row,-(1:4)] <- 
            c(test_summ[,c("Mean", "SD")], 
              p = ifelse(bayesp > .5, 1 - bayesp, bayesp) * 2,  # two-sided Bayesian p values
              test_summ[,2 + seq_along(CI_SIZE)])
        }
      }
    }
  }
  
  # don't save T2 if complete pooling was used ----
  # Why? I think it would be worthwhile
  if (pooling != "complete"){
    result_row$gof[[1]] <- tibble::add_row(result_row$gof[[1]])   # T1 & T2
    result_row$gof[[1]][2,-(1:2)] <- aggregate_ppp(gof_group, stat = "T2")
  }
  result_row$gof[[1]]$type <- c("T1", if(pooling!="complete"){"T2"})
  result_row$gof[[1]]$focus <- c("mean", if(pooling!="complete"){"cov"})
  
  result_row$gof[[1]][1,-(1:2)] <- aggregate_ppp(gof_group)
  
  estimation_time <- unlist(estimation_time)
  result_row$estimation[[1]] <- tibble::tibble(
    condition = names(estimation_time)
    , time_difference = unname(estimation_time)
  )

  
  # save model objects to the working directory if requested by user ----
  if(all_options$save_models){
    save(treebugs_fit, file = paste0(
      paste(
        c(
          gsub(model, pattern = ".eqn|.EQN", replacement = "")
          , gsub(dataset, pattern = ".csv|.CSV", replacement = "")
          , pooling
          , method
        )
        , collapse = "_"
      )
      , ".RData"
      )
    )
  }
  # return ----
  result_row
}



mpt_treebugs_safe <- purrr::possibly(
  .f = mpt_treebugs
  , otherwise = list()
  , quiet = FALSE
)
