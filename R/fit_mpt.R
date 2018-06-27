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
  
  
  
  # extraneous columns, etc. should be removed here
  
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
  
    running_on_windows <- Sys.info()[["sysname"]]=="Windows"
  
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