# plot_results <- function (results, save = TRUE, write.csv = TRUE){
#   
#   shapes <- c(16, 18, 15, 1, 0, 8, 11, 12)
# 
#   prefix <- paste0(gsub("\\.eqn", "", results$model[1]), "_", 
#                             gsub("\\.", "_", paste0(results$dataset[1],"_")))
#   
#   if (write.csv){
#     readr::write_csv(tidyr::unnest(results, est_group), paste0(prefix,"estimates.csv"))
#     readr::write_csv(tidyr::unnest(results, gof), paste0(prefix,"gof.csv"))
#   }
#   
#   
#   dd <- ggplot2::position_dodge(w = .75)
#   gg_est1 <- tidyr::unnest(results, est_group) %>%
#     ggplot2::ggplot(ggplot2::aes(y = est, x = parameter, 
#                col=interaction(method, pooling,package),
#                shape=interaction(method, pooling, package))) +
#     ggplot2::facet_grid(.~condition) +
#     ggplot2::geom_errorbar(ggplot2::aes(ymin = est-se, ymax = est+se), position = dd, 
#                   width = 0.4)+
#     ggplot2::geom_point(position = dd) + ylim(0,1) + 
#     ggplot2::scale_shape_manual(values=shapes) +
#     ggplot2::theme_bw()
#   plot(gg_est1)
#   if(save) ggplot2::ggsave(paste0(prefix,"estimates.pdf"), gg_est1, h = 4.5, w = 10)
#   
#   
#   res_between <-  tidyr::unnest(results, test_between) 
#   if (nrow(res_between) > 0){
#     
#     if (write.csv) readr::write_csv(tidyr::unnest(results, test_between), paste0(prefix,"test_between.csv"))
#     
#     gg_est2 <- ggplot2::ggplot(
#       res_between
#       , ggplot2::aes(
#         y = est_diff
#         , x = parameter
#         , col = interaction(method, pooling, package)
#         , shape = interaction(method, pooling, package)
#       )
#     ) +
#     ggplot2::facet_grid(condition2~condition1) +
#     ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_0.025, ymax = ci_0.975), 
#                   position = dd, width = 0.5) +
#     ggplot2::geom_point(position = dd) + ylim(-1,1) + 
#     ggplot2::scale_shape_manual(values=shapes) +
#     ggplot2::theme_bw() + 
#     ggplot2::geom_hline(yintercept = 0, lty = 2)
#     
#     if(save) ggplot2::ggsave(paste0(prefix,"test_between.pdf"), gg_est2, h = 4.5, w = 8)
#   }
#   
#   
#   gg_gof1 <-  tidyr::unnest(results, gof) %>%
#     # filter(focus == "mean") %>%
#     ggplot2::ggplot(ggplot2::aes(y = p, 
#                x = interaction(method, pooling, package))) + 
#     ggplot2::geom_point() + ylim(0, 1) + 
#     ggplot2::geom_hline(yintercept = .05, lty = 2)+
#     ggplot2::theme_bw() + ggplot2::coord_flip() +
#     ggplot2::facet_wrap(~focus) +
#     ggplot2::ggtitle("Goodness of fit")
#   plot(gg_gof1)
#   if(save) ggplot2::ggsave(paste0(prefix,"gof.pdf"), gg_gof1, h = 4, w = 6)
#   
#   
#   if (nrow(res_between) > 0){
#     if (write.csv) readr::write_csv(tidyr::unnest(results, gof_group), paste0(prefix,"gof_group.csv"))
#     
#     gg_gof2 <- tidyr::unnest(results, gof_group) %>%
#       ggplot(ggplot2::aes(y = p, 
#                  x = interaction(method, pooling, package), 
#                  col = condition)) + 
#       ggplot2::geom_point() + ylim(0, 1) + 
#       ggplot2::geom_hline(yintercept = .05, lty = 2)+
#       ggplot2::theme_bw() +
#       ggplot2::coord_flip() +
#       ggplot2::facet_wrap(~focus) +
#       ggplot2::ggtitle("Goodness of fit")
#     plot(gg_gof2)
#     if(save) ggplot2::ggsave(paste0(prefix,"gof_group.pdf"), gg_gof2, h = 4, w = 8)
#   }
# }


#' Plot multiverseMPT
#' 
#' Plot the results from a multiverse MPT analysis.
#'  
#' @param x An object of class \code{multiverseMPT}.
#' @param which Character. Which information should be plotted? Possible
#' values are
#' \code{"est"} for parameter estimates,
#' \code{"test_between"} for between-subjects comparisions,
#' \code{"gof1"} for overall goodness-of-fit statistics, and
#' \code{"gof2"} for group-wise goodness-of-fit statistics.
#' @param save Logical. Indicates whether the plot should also be saved as a .pdf file.
#' @param ... ignored.
#' 
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom graphics plot
#' @export

plot.multiverseMPT <- function(x, which = "est", save = FALSE, ...){
  
  shapes <- c(16, 18, 15, 1, 0, 8, 11, 12, 4, 6)
  
  results <- x
  prefix <- paste0(gsub("\\.eqn", "", results$model[1]), "_", 
                   gsub("\\.", "_", paste0(results$dataset[1],"_")))
  
  # if (write.csv){
  #   readr::write_csv(tidyr::unnest(results, .data$est_group), paste0(prefix,"estimates.csv"))
  #   readr::write_csv(tidyr::unnest(results, .data$gof), paste0(prefix,"gof.csv"))
  # }
  
  
  dd <- ggplot2::position_dodge(width = .75)
  
  est_group <- tidyr::unnest(data = results, .data$est_group)
  est_group$approach <- interaction(est_group$method, est_group$pooling, est_group$package)
  
  gg_est1 <-
    ggplot2::ggplot(est_group) +
    ggplot2::aes_(
      y = ~ est
      , x = ~ parameter
      , col = ~ approach
      # , shape = shapes
    ) +
    ggplot2::facet_grid(facets = ". ~ condition") +
    ggplot2::geom_errorbar(ggplot2::aes_(ymin = ~ci_0.025, ymax = ~ci_0.975), position = dd, width = 0.4) +
    ggplot2::geom_point(position = dd) + 
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    # ggplot2::ylim(0, 1) + 
    ggplot2::scale_shape_manual(values = shapes)

  if(save) ggplot2::ggsave(paste0(prefix, "estimates.pdf"), gg_est1, h = 4.5, w = 10)
  
  if("est" %in% which)
    return(gg_est1)
  
  res_between <-  tidyr::unnest(results, .data$test_between)
  res_between$approach <- interaction(res_between$method, res_between$pooling, res_between$package)
  
  if (nrow(res_between) > 0){
    # if (write.csv) readr::write_csv(res_between, paste0(prefix,"test_between.csv"))
    
    
    gg_est2 <- ggplot2::ggplot(
      res_between
      , ggplot2::aes_(
        y = ~ est_diff, x = ~ parameter
        , col = ~ approach
        , shape = ~ approach
      )
    ) +
    ggplot2::facet_grid("condition2 ~ condition1") +
    ggplot2::geom_errorbar(
      ggplot2::aes_(
        ymin = ~ ci_0.025
        , ymax = ~ ci_0.975
      )
      , position = dd, width = 0.5
    ) +
    ggplot2::geom_point(position = dd) + ggplot2::ylim(-1, 1) + 
    ggplot2::scale_shape_manual(values=shapes) +
    ggplot2::geom_hline(yintercept = 0, lty = 2)

    if(save) ggplot2::ggsave(paste0(prefix,"test_between.pdf"), gg_est2, h = 4.5, w = 8)
    
    if("test_between" %in% which){
      return(gg_est2)
    }
  }
  
  gof <- tidyr::unnest(results, .data$gof)
  gof$approach <- interaction(gof$method, gof$pooling, gof$package)
  
  gg_gof1 <-
    # filter(focus == "mean") %>%
    ggplot2::ggplot(
      gof, 
      ggplot2::aes_(y = ~ p, x = ~ approach)
    ) + 
    ggplot2::geom_point() + 
    ggplot2::ylim(0, 1) + 
    ggplot2::geom_hline(yintercept = .05, lty = 2)+
    ggplot2::coord_flip() +
    ggplot2::facet_wrap( ~ focus) +
    ggplot2::ggtitle("Goodness of fit")
  
  if(save) ggplot2::ggsave(paste0(prefix,"gof.pdf"), gg_gof1, h = 4, w = 6)
  
  if("gof1" %in% which)
    return(gg_gof1)
  
  
  if (nrow(res_between) > 0){
    # if (write.csv) readr::write_csv(tidyr::unnest(results, .data$gof_group), paste0(prefix,"gof_group.csv"))
    
    gof_group <- tidyr::unnest(results, .data$gof_group)
    gof_group$approach <- interaction(gof_group$method, gof_group$pooling, gof_group$package)
    
    gg_gof2 <-
      ggplot2::ggplot(
        gof_group
        , ggplot2::aes_(
          y = ~ p
          , x = ~ approach
          , col = ~ condition)
        ) + 
      ggplot2::geom_point() + ggplot2::ylim(0, 1) + 
      ggplot2::geom_hline(yintercept = .05, lty = 2)+
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~ focus) +
      ggplot2::ggtitle("Goodness of fit")
    if(save) ggplot2::ggsave(paste0(prefix, "gof_group.pdf"), gg_gof2, h = 4, w = 8)
    
    if("gof2" %in% which)
      return(gg_gof2)
  }
}



# x <- results
# 
# 
# x <- tidyr::tidyr::unnest(x, est_group)
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
