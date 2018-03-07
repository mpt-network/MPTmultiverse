#' @export

plot_results <- function (results, save = TRUE, write.csv = TRUE){
  
  shapes <- c(16, 18, 15, 1, 0, 8, 11, 12)

  prefix <- paste0(gsub("\\.eqn", "", results$model[1]), "_", 
                            gsub("\\.", "_", paste0(results$dataset[1],"_")))
  
  if (write.csv){
    write_csv(unnest(results, est_group), paste0(prefix,"estimates.csv"))
    write_csv(unnest(results, gof), paste0(prefix,"gof.csv"))
  }
  
  
  dd <- position_dodge(w = .75)
  gg_est1 <- unnest(results, est_group) %>%
    ggplot(aes(y = est, x = parameter, 
               col=interaction(method, pooling,package),
               shape=interaction(method, pooling,package))) +
    facet_grid(.~condition) +
    geom_errorbar(aes(ymin = est-se, ymax = est+se), position = dd, 
                  width = 0.4)+
    geom_point(position = dd) + ylim(0,1) + 
    scale_shape_manual(values=shapes) +
    theme_bw()
  plot(gg_est1)
  if(save) ggsave(paste0(prefix,"estimates.pdf"), gg_est1, h = 4.5, w = 10)
  
  
  res_between <-  unnest(results, test_between) 
  if (nrow(res_between) > 0){
    if (write.csv) write_csv(unnest(results, test_between), paste0(prefix,"test_between.csv"))
    
    gg_est2 <- ggplot(res_between, aes(y = est_diff, x = parameter, 
                                       col=interaction(method, pooling,package),
                                       shape=interaction(method, pooling,package))) +
      facet_grid(condition2~condition1) +
      geom_errorbar(aes(ymin = ci_0.025, ymax = ci_0.975), 
                    position = dd, width = 0.5)+
      geom_point(position = dd) + ylim(-1,1) + 
      scale_shape_manual(values=shapes) +
      theme_bw() + geom_hline(yintercept = 0, lty = 2)
    plot(gg_est2)
    if(save) ggsave(paste0(prefix,"test_between.pdf"), gg_est2, h = 4.5, w = 8)
  }
  
  
  gg_gof1 <-  unnest(results, gof) %>%
    # filter(focus == "mean") %>%
    ggplot(aes(y = p, 
               x = interaction(method, pooling, package))) + 
    geom_point() + ylim(0, 1) + 
    geom_hline(yintercept = .05, lty = 2)+
    theme_bw() + coord_flip() +
    facet_wrap(~focus) +
    ggtitle("Goodness of fit")
  plot(gg_gof1)
  if(save) ggsave(paste0(prefix,"gof.pdf"), gg_gof1, h = 4, w = 6)
  
  
  if (nrow(res_between) > 0){
    if (write.csv) write_csv(unnest(results, gof_group), paste0(prefix,"gof_group.csv"))
    
    gg_gof2 <- unnest(results, gof_group) %>%
      ggplot(aes(y = p, 
                 x = interaction(method, pooling, package), 
                 col = condition)) + 
      geom_point() + ylim(0, 1) + 
      geom_hline(yintercept = .05, lty = 2)+
      theme_bw() +
      coord_flip() +
      facet_wrap(~focus) +
      ggtitle("Goodness of fit")
    plot(gg_gof2)
    if(save) ggsave(paste0(prefix,"gof_group.pdf"), gg_gof2, h = 4, w = 8)
  }
}

#' @export

plot.multiverseMPT <- function(x, which = "est", save = FALSE, write.csv = FALSE){
  
  if("est" %in% which){
    
  }

  shapes <- c(16, 18, 15, 1, 0, 8, 11, 12)
  
  prefix <- paste0(gsub("\\.eqn", "", results$model[1]), "_", 
                   gsub("\\.", "_", paste0(results$dataset[1],"_")))
  
  if (write.csv){
    write_csv(unnest(results, est_group), paste0(prefix,"estimates.csv"))
    write_csv(unnest(results, gof), paste0(prefix,"gof.csv"))
  }
  
  
  dd <- position_dodge(w = .75)
  
  
  gg_est1 <- unnest(results, est_group) %>%
    ggplot(aes(y = est, x = parameter, 
               col=interaction(method, pooling,package),
               shape=interaction(method, pooling, package))) +
    facet_grid(.~condition) +
    geom_errorbar(aes(ymin = est-se, ymax = est+se), position = dd, 
                  width = 0.4)+
    geom_point(position = dd) + ylim(0,1) + 
    scale_shape_manual(values=shapes)

  if(save) ggsave(paste0(prefix,"estimates.pdf"), gg_est1, h = 4.5, w = 10)
  
  if("est" %in% which)
    return(gg_est1)
  
  res_between <-  unnest(results, test_between)
  
  
  if (nrow(res_between) > 0){
    if (write.csv) write_csv(unnest(results, test_between), paste0(prefix,"test_between.csv"))
    
    
    gg_est2 <- ggplot(res_between, aes(y = est_diff, x = parameter, 
                                       col=interaction(method, pooling,package),
                                       shape=interaction(method, pooling,package))) +
      facet_grid(condition2~condition1) +
      geom_errorbar(aes(ymin = ci_0.025, ymax = ci_0.975), 
                    position = dd, width = 0.5)+
      geom_point(position = dd) + ylim(-1,1) + 
      scale_shape_manual(values=shapes) +
      geom_hline(yintercept = 0, lty = 2)

    if(save) ggsave(paste0(prefix,"test_between.pdf"), gg_est2, h = 4.5, w = 8)
    
    if("test_between" %in% which){
      return(gg_est2)
    }
  }
  
  
  gg_gof1 <-  unnest(results, gof) %>%
    # filter(focus == "mean") %>%
    ggplot(aes(y = p, 
               x = interaction(method, pooling, package))) + 
    geom_point() + ylim(0, 1) + 
    geom_hline(yintercept = .05, lty = 2)+
    coord_flip() +
    facet_wrap(~focus) +
    ggtitle("Goodness of fit")
  
  if(save) ggsave(paste0(prefix,"gof.pdf"), gg_gof1, h = 4, w = 6)
  
  if("gof1" %in% which)
    return(gg_gof1)
  
  
  if (nrow(res_between) > 0){
    if (write.csv) write_csv(unnest(results, gof_group), paste0(prefix,"gof_group.csv"))
    
    gg_gof2 <- unnest(results, gof_group) %>%
      ggplot(aes(y = p, 
                 x = interaction(method, pooling, package), 
                 col = condition)) + 
      geom_point() + ylim(0, 1) + 
      geom_hline(yintercept = .05, lty = 2)+
      coord_flip() +
      facet_wrap(~focus) +
      ggtitle("Goodness of fit")
    if(save) ggsave(paste0(prefix,"gof_group.pdf"), gg_gof2, h = 4, w = 8)
    
    if("gof2" %in% which)
      return(gg_gof2)
  }
}



# x <- results
# 
# 
# x <- tidyr::unnest(x, est_group)
# 
# y.values <- data.frame(
#   condition = x$condition
#   , method = paste0(x$method, "_", x$pooling, "_pooling")
#   , parameter = x$parameter
#   , tendency = x$est
#   , dispersion = x$se
#   , lower_limit = x$ci_0.025
#   , upper_limit = x$ci_0.025
# )
# 
# # re-order factor
# y.values$method <- factor(
#   y.values$method
#   , levels = c(
#     "asymptotic_complete_pooling"
#     , "simple_complete_pooling"
#     , "asymptotic_no_pooling"
#     , "PB/MLE_no_pooling"
#     , "simple_no_pooling"
#     , "beta_partial_pooling"
#     , "trait_uncorrelated_partial_pooling"
#     , "trait_partial_pooling"
#   )
# )
# 
# x$method <- paste0(x$method, "_", x$pooling, "_pooling")
# 
# par(mar = c(5, 4, 4, 14), las = 1)
# plot_options <- papaja:::apa_factorial_plot_single(
#   aggregated = x
#   , y.values = y.values
#   , id = "id"
#   , dv = "est"
#   , factors = c("parameter", "method", "condition")
#   , ylim = c(0, 1)
#   , plot = c("points", "error_bars")
#   , args_points = list(
#     cex = .8
#     , pch = c(rep(21, 2), rep(22, 3), rep(23, 3))
#     , bg = c("seagreen4", "indianred3", "skyblue")
#   )
#   , args_error_bars = list(length = .01)
#   , jit = .5
#   , args_legend = list(x = "right", inset = -.5, xpd = NA)
# )
