#' Check results from a multiverse analysis
#' 
#' This is a helper function to see if the model estimation worked as intended.
#' 
#' @param results An object of class multiverseMPT.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export

check_results <- function(results) {
  #browser()
  expected <- structure(list(
    pooling = c("no", "no", "no", "complete", "no", "complete", "partial", 
                "partial", "partial", "partial"), 
    package = c("MPTinR", "MPTinR", "MPTinR", "MPTinR", "TreeBUGS", "TreeBUGS", 
                "TreeBUGS", "TreeBUGS", "TreeBUGS", "TreeBUGS"), 
    method = c("NPB/MLE", "PB/MLE", "asymptotic", "asymptotic", "simple", 
               "simple", "trait", "trait_uncorrelated", "beta", "betacpp")), 
    .Names = c("pooling", "package", "method"), 
    class = c("tbl_df", "tbl", "data.frame"
    ), row.names = c(NA, -10L))
  missing <- dplyr::anti_join(expected, results[, 3:5], by = c("pooling", "package", "method"))
  if (nrow(missing) > 0) {
    cat("## Following analysis approaches missing from results:\n", 
        paste(apply(missing, 1, paste, collapse = ", "), collapse = "\n"), 
        "\n\n\n")
  }
  
  ### MPTinR: no pooling ###
  
  cat("## MPTinR: no pooling\n")
  
  tryCatch({
    for(meth in c("asymptotic", "PB/MLE", "NPB/MLE")){
      
      conv_mptinr_no <- results %>% 
        dplyr::filter(.data$package == "MPTinR" & .data$pooling == "no" & .data$method == meth) %>% 
        dplyr::select("convergence") %>% 
        tidyr::unnest()
      
      not_id <- results %>% 
        dplyr::filter(.data$package == "MPTinR" & .data$pooling == "no" & .data$method == meth) %>% 
        dplyr::select("est_indiv") %>% 
        tidyr::unnest() %>% 
        dplyr::group_by(.data$condition) %>% 
        dplyr::summarise(proportion = mean(!.data$identifiable))
      
      not_id2 <- results %>% 
        dplyr::filter(.data$package == "MPTinR" & .data$pooling == "no" & .data$method == meth) %>% 
        dplyr::select("est_indiv") %>% 
        tidyr::unnest() %>%
        dplyr::filter(!.data$identifiable) %>% 
        dplyr::group_by(.data$condition) %>% 
        dplyr::summarise(not_identified = list(broom::tidy(table(.data$parameter)))) %>% 
        tidyr::unnest(.data$not_identified) %>% 
        suppressWarnings()
      
      if (any(not_id$proportion > 0)) {
        cat("Based on", meth, "method, proportion of participants with non-identified parameters:\n")
        cat(format(not_id)[-c(1,3)], "", sep = "\n")
        
        cat("Based on", meth, "CIs, table of non-identified parameters:\n")
        cat(format(not_id2)[-c(1,3)], sep = "\n")
        
      } else {
        cat("Based on", meth, "CIs, all parameters of all participants seem to be identifiable.\n")
      }
      cat("\n")
    }
  }, error = function(e) 
    cat("Convergence checks failed for unkown reason.\n"))
  
  cat("\n")
  
  
  ### MPTinR: complete pooling ###
  
  cat("## MPTinR: complete pooling\n")
  
  tryCatch({
    conv_mptinr_comp <- results %>%
      dplyr::filter(.data$package == "MPTinR" & .data$pooling == "complete") %>%
      dplyr::select("convergence") %>%
      tidyr::unnest()
    
    comp_prob <- (conv_mptinr_comp$convergence != 0) | 
      (conv_mptinr_comp$rank.fisher != conv_mptinr_comp$n.parameters)
    
    if (any(comp_prob)) {
      cat("Convergence problems:\n")
      cat(format(conv_mptinr_comp[comp_prob,])[-c(1,3)], "", sep = "\n")
    } else {
      cat("No convergence problems.\n")
    }
  }, error = function(e) 
    cat("Convergence checks failed for unkown reason.\n"))
  
  cat("\n\n")
  
  ### TreeBUGS
  res_tree <- results %>% 
    dplyr::filter(.data$package == "TreeBUGS") %>% 
    dplyr::select(!!c("pooling", "package", "method", "convergence"))
  
  for (i in seq_len(nrow(res_tree))) {
    cat("## ", paste(res_tree[i, c(2,1,3)], collapse = ", "), ":\n", sep = "")
    
    tmp_convergence <- res_tree[i, ]$convergence[[1]] %>% 
      dplyr::filter(.data$Rhat > getOption("MPTmultiverse")$treebugs$Rhat_max) 
    
    if(nrow(tmp_convergence) > 0) {
      cat(nrow(tmp_convergence), "parameters with Rhat >", getOption("MPTmultiverse")$treebugs$Rhat_max, ":\n")
      cat(paste(tmp_convergence$parameter, collapse = ", "))
    } else {
      cat("All Rhat <", getOption("MPTmultiverse")$treebugs$Rhat_max, ".\n")
    }
    
    tmp_neff <- res_tree[i,]$convergence[[1]] %>% 
      dplyr::filter(!is.na(.data$Rhat), .data$n.eff < getOption("MPTmultiverse")$treebugs$Neff_min)
    
    if(nrow(tmp_neff) > 0) {
      cat(nrow(tmp_neff), "parameters with effect sample size n.eff <", 
          getOption("MPTmultiverse")$treebugs$Neff_min, ":\n")
      cat(paste(tmp_neff$parameter, collapse = ", "))
    } else {
      cat("All effect sample sizes >", getOption("MPTmultiverse")$treebugs$Neff_min, ".\n")
    }
    
    cat("\n\n")
  }
  
  
}

#' Write check_results
#'
#' Helper function to write the output from check_results() to a file.
#'
#' @param DATA_FILE Character. File name to use.
#' @param results An object of class multiverseMPT.
#' @export

write_check_results <- function(DATA_FILE, results){
  sink(paste0(DATA_FILE, "_check_results.txt"))
  cat("################ OPTIONS ################\n\n")
  cat("TreeBUGS:\n") ; print(getOption("MPTmultiverse")$treebugs)
  cat("\nMPTinR:\n") ; print(getOption("MPTmultiverse")$mptinr)
  cat("\nCI_SIZE: ", getOption("MPTmultiverse")$ci_size, "\n")
  cat("MAX_CI_INDIV = ", getOption("MPTmultiverse")$max_ci_indiv, "\n\n")
  cat("################ CHECK RESULTS ################\n\n")
  print(check_results(results))
  sink()
}
