#' @export

compare <- function(x, ...) {
  UseMethod("compare")
}

#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom graphics plot
#' @export

compare.multiverseMPT <- function(
  x
  , ...
) {
  
  args <- list(...)
  
  if(is.null(args$which)) {
    args$which <- "est_group"
  }
  
  results <- tidyr::unnest(results, .data[[args$which]])
  
  results$inter <- interaction(
    results$method
    , results$pooling
    , results$package
    , drop = TRUE
    , sep = " "
  )
  
  results <- results %>%
    dplyr::group_by(.data$model, .data$dataset, .data$inter, .data$condition, .data$parameter) %>%
    dplyr::mutate(within = seq_along(.data$parameter)) %>%
    dplyr::ungroup()
  
  # Calculate CCC of all pairwise combinations
  pairs <- utils::combn(sort(levels(results$inter)), 2)
  all_pars_list <- vector("list", ncol(pairs))
  
  for (i in seq_len(ncol(pairs))) {
    tmp_dat <- results %>% 
      dplyr::filter(.data$inter %in% pairs[, i]) %>% 
      dplyr::select(.data$model, .data$dataset, .data$condition, .data$within, .data$parameter, .data$est, .data$inter) %>% 
      dplyr::group_by(.data$model, .data$dataset) %>% 
      dplyr::mutate(n_conditions = length(unique(.data$condition))) %>% 
      dplyr::ungroup() %>%
      tidyr::spread(key = .data$inter, value = .data$est)
    colnames(tmp_dat)[(ncol(tmp_dat)-1):ncol(tmp_dat)] <- c("x", "y")
    tmp_dat$cond_x <- pairs[1, i]
    tmp_dat$cond_y <- pairs[2, i]
    all_pars_list[[i]] <- tmp_dat 
  }
  
  all_pars <- dplyr::bind_rows(all_pars_list)
  all_pars
}
