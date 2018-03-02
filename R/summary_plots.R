library("ggplot2")

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

