
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

#' @export

mpt_treebugs <- function (method, dataset, data, model,
                          col_id = "id", col_condition = "condition"){
  all_options <- getOption("mpt.comparison")
  
  TREEBUGS_MCMC <- all_options$treebugs
  CI_SIZE <- all_options$ci_size
  
  # dlist <- prepare_data(model, data, col_id = "id", col_condition = "condition")
  conditions <- levels(factor(data[[col_condition]]))
  parameters <- check.mpt(model)$parameters
  col_freq <- get_eqn_categories(model)

  data$id <- data[, col_id]
  data$condition <- data[,col_condition]
  freq_list <- split(data[, col_freq], f = data[, col_condition])
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
                                 parameters = parameters)
  if (method == "simple_pooling"){
    method <- "simple"
    
    # pooling: aggregate across participants
    data <- aggregate(data[,col_freq], list(condition = data$condition), sum)
    data[[col_condition]] <- data$condition
    data[[col_id]] <- data$id <- 1:nrow(data)
    if(col_condition!="condition"){
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
  for (i in seq_along(conditions)){
    cond <- factor(conditions[i], conditions)
    sel_condition <- data[[col_condition]] == conditions[i]
    data_group <- data[sel_condition, col_freq]   #freq_list[[i]]
    rownames(data_group) <- data[[col_id]][sel_condition]
    
    fit_args <- list(eqnfile=model,
                     data = data_group,
                     n.chain = TREEBUGS_MCMC$n.chain,
                     n.iter = TREEBUGS_MCMC$n.iter,
                     n.adapt = TREEBUGS_MCMC$n.adapt,
                     n.burnin = TREEBUGS_MCMC$n.burnin,
                     n.thin = TREEBUGS_MCMC$n.thin)
    if (method == "simple"){
      fit_args["n.adapt"] <- NULL
      fit_args <- c(fit_args, cores = unname(TREEBUGS_MCMC$n.CPU))
    }
    # print(c(fit_args, prior_args))
    treebugs_fit[[i]] <- do.call(paste0(method, "MPT"), args = c(fit_args, prior_args))
    summ <- treebugs_fit[[i]]$mcmc.summ
    
    # continue MCMC sampling (only for betaMPT and traitMPT)
    ext_cnt <- 0
    try({
      while (
        ext_cnt < TREEBUGS_MCMC$extend_max && method %in% c("beta", "trait") &&
        (any(na.omit(summ[,"Rhat"]) > TREEBUGS_MCMC$Rhat_max)  ||
         any(summ[summ[,"n.eff"] > 0,"n.eff"] < TREEBUGS_MCMC$Neff_min, na.rm = TRUE)) ){
        cat("Drawing additional samples for method = ", method, 
            ". max(Rhat) = ", round(max(na.omit(summ[summ[,"Rhat"] > 0,"Rhat"])), 2),
            " ; min(n.eff) = ", round(min(summ[summ[,"n.eff"] > 0,"n.eff"], na.rm = TRUE), 1), "\n")
        
        treebugs_fit[[i]] <- extendMPT(treebugs_fit[[i]],
                                       n.iter = TREEBUGS_MCMC$n.iter,
                                       n.thin = TREEBUGS_MCMC$n.thin,
                                       n.adapt = TREEBUGS_MCMC$n.adapt)
        summ <- treebugs_fit[[i]]$mcmc.summ
        ext_cnt <- ext_cnt + 1
      }
    })
    
    # convergence summary (n.eff / Rhat / all estimates)
    tsum <- as_tibble(summ) %>% 
      mutate(parameter = rownames(summ),
             condition = as.character(cond)) %>% 
      select(condition, parameter, Mean:Rhat)
    result_row$convergence[[1]] <- bind_rows(result_row$convergence[[1]], tsum)
    
    # parameter estimates
    summMPT <- summarizeMPT(treebugs_fit[[i]]$runjags$mcmc,
                            mptInfo = treebugs_fit[[i]]$mptInfo,
                            probs = CI_SIZE)
    
    sel_group <- result_row$est_group[[1]]$condition == conditions[i]
    result_row$est_group[[1]][sel_group,-(1:2)] <-
      summMPT$groupParameters$mean[paste0("mean_", parameters),1:6]
    
    if (pooling != "complete"){
      # # old: array filled into data frame
      # result_row$est_indiv[[1]][sel_ind,-(1:3)] <-
      #   summMPT$individParameters[parameters,,1:(2+length(CI_SIZE))]
      sel_ind <- result_row$est_indiv[[1]]$condition == conditions[i]
      dimnames(summMPT$individParameters)$ID <- rownames(data_group)
      tmp <- summMPT$individParameters[parameters,,1:(2+length(CI_SIZE)), drop = FALSE] %>%
        melt %>% 
        spread("Statistic", "value")
      colnames(tmp) <- c("parameter", "id", colnames(result_row$est_indiv[[1]])[-(1:3)])
      tmp[[col_condition]] <- cond
      result_row$est_indiv[[1]][sel_ind,] <-
        left_join(result_row$est_indiv[[1]][sel_ind,] %>%
                    select("id", "condition", "parameter"),
                  tmp, by = c("parameter", "id", condition = col_condition))
    }
    
    gof_group[[i]] <- PPP(treebugs_fit[[i]], M = TREEBUGS_MCMC$n.PPP, type = "G2",
                          T2 = pooling != "complete", nCPU = TREEBUGS_MCMC$n.CPU)
    
    sel_gof <- result_row$gof_group[[1]]$condition == conditions[i]
    result_row$gof_group[[1]][sel_gof,] <-
      result_row$gof_group[[1]] %>%
      filter(condition == conditions[i]) %>%
      mutate(condition = conditions[i],
             type = "T1_G2", focus = "mean",
             stat_obs = mean(gof_group[[i]]$T1.obs),
             stat_pred = mean(gof_group[[i]]$T1.pred),
             p = gof_group[[i]]$T1.p)
    
    if (pooling != "complete"){
      result_row$gof_group[[1]] <- add_row(result_row$gof_group[[1]],
                                           condition = cond,
                                           type = "T2", focus = "cov",
                                           stat_obs = mean(gof_group[[i]]$T2.obs),
                                           stat_pred = mean(gof_group[[i]]$T2.pred),
                                           p = gof_group[[i]]$T2.p)
      
      sel_fog_ind <- result_row$gof_indiv[[1]]$condition == conditions[i]
      result_row$gof_indiv[[1]][sel_fog_ind,] <-
        result_row$gof_indiv[[1]][sel_fog_ind,] %>%
        mutate(condition = conditions[i],
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
          test_between <- betweenSubjectMPT(treebugs_fit[[i]], treebugs_fit[[j]], 
                                            par1 = p, stat = "x-y")
          test_summ <- summarizeMCMC(test_between$mcmc, probs = CI_SIZE)
          bayesp <- mean(do.call("rbind", test_between$mcmc) <= 0)
          
          sel_row <- 
            result_row$test_between[[1]]$parameter == p &
            result_row$test_between[[1]]$condition1 == conditions[i] &
            result_row$test_between[[1]]$condition2 == conditions[j]
          
          result_row$test_between[[1]][sel_row,-(1:3)] <- 
            c(test_summ[,c("Mean", "SD")], 
              p = ifelse(bayesp > .5, 1 - bayesp, bayesp) * 2,  # two-sided Bayesian p values
              test_summ[,2 + seq_along(CI_SIZE)])
        }
      }
    }
  }
  
  result_row$gof[[1]] <- add_row(result_row$gof[[1]])   # T1 & T2
  result_row$gof[[1]]$type <- c("T1_G2", "T2")
  result_row$gof[[1]]$focus <- c("mean", "cov")
  result_row$gof[[1]][1,-(1:2)] <- aggregate_ppp(gof_group)
  if (pooling != "complete")
    result_row$gof[[1]][2,-(1:2)] <- aggregate_ppp(gof_group, stat = "T2")
  
  if(all_options$save_models){
    save(treebugs_fit, file = paste0(model, data, method, ".RData"))
  }
  result_row
}


#' @export

mpt_treebugs_safe <- possibly(mpt_treebugs, otherwise = list())
