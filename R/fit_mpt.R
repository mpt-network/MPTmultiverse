#' Main function
#'
#' Does a lot of nice stuff
#'
#' @param method Character. A vector specifying which analysis approaches should be performed.
#' @param dataset Character.
#' @param data A \code{data.frame} containing the data.
#' @param model A model definition, typically the path to an \code{.eqn} file.
#' @param id Character. Name of the column that contains the subject identifier.
#'   If not specified, it is assumed that each row represents observations from one participant.
#' @param condition Character. Name of the column specifying a between-subjects factor.
#'   If not specified, no between-subjects comparisons are performed.
#'
#' @export

fit_mpt <- function(
  method
  , dataset
  , data
  , model
  , id = NULL
  , condition = NULL
) {
  
  
  # set options ----
  silent <- getOption("MPTmultiverse")$silent
  
  runjags::runjags.options(
    silent.jags = silent
    , silent.runjags = silent
  )
  
  # prepare data ----
  if(is.null(condition)) {
    data$ExpCond <- factor("no_condition")
    condition <- "ExpCond"
  }
  
  if(is.null(id)) {
    data$Subject <- 1:nrow(data)
    id <- "Subject"
  }
  
  # remove extraneous colums ----
  valid_cols <- c(id, condition, get_eqn_categories(model))
  data <- data[, valid_cols]
  
  # Ensure that id and condition are factors, also drops unused levels
  # data[[id]] <- factor(data[[id]])
  data[[condition]] <- factor(data[[condition]])

  
  res <- list()
  
  # MPTinR part ----
  res[["mptinr"]] <- mpt_mptinr(
    dataset = dataset
    , data = data
    , model = model
    , method = intersect(method, c("asymptotic_complete", "asymptotic_no", "pb_no"))
    , id = id
    , condition = condition
  )
  
  # HMMTreeR part ----
  if("latent_class" %in% method) {
    
    # HMMTreeR installed?
    if(suppressWarnings(requireNamespace("HMMTreeR"))) {
  
    running_on_windows <- Sys.info()[["sysname"]]=="Windows"
    
      # ensure that no fixed parameter values are present
      test <- 1
      test <- tryCatch(simplify_eqn(
        model_filename = model
        , eqn_filename = "tmp_eqn_HMMTree.eqn"
        , data = data
        , id = id
        , condition = condition
      ))
      file.remove("tmp_eqn_HMMTree.eqn")

      if(test==0) {
    
        if(running_on_windows) {
          res[["hmmtreer"]] <- fit_lc(
            dataset = dataset
            , data = data
            , model = model
            , id = id
            , condition = condition
          )
        } else {
          message("Latent-class multinomial models can currently only be estimated on Windows -- sorry.")
        }
      } else {
        message("The specified .eqn file seems to contain fixed parameter values:
                The current implementation of latent-class models does not support this type of .eqn files.
                Therefore, latant-class models are not estimated.")
      }
    } else {
      message("HMMTreeR is not installed on your system. Therefore, latent-class analyses are skipped.")
    }
  }
  
  # TreeBUGS part ----
  res[["treebugs"]] <- dplyr::bind_rows(
    purrr::map(
      intersect(method, c("simple", "simple_pooling", "trait", "beta", "trait_uncorrelated"))
      , mpt_treebugs
      , dataset = dataset
      , data = data
      , model = model
      , id = id
      , condition = condition
    )
  )

  y <- dplyr::bind_rows(res)
  class(y) <- c("multiverseMPT", class(y))
  y
}