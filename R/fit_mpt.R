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
  
  # check MPT file
  mpt_model <- TreeBUGS::readEQN(model)
  
  if(!is.data.frame(mpt_model)) {
    "I can't comprehend your .eqn file."
  }

  
  # remove extraneous colums and check if all specified columns are present
  # in data
  freq_cols <- get_eqn_categories(model)
  valid_cols <- c(id, condition, freq_cols)
  check_cols <- valid_cols %in% colnames(data)
  
  if(!all(check_cols)) {
    stop("Variable \"", paste(valid_cols[!check_cols], collapse = ", "), "\" not found in data.frame.")
  }
  
  data <- data[, valid_cols]
  
  
  # Check NAs ----
  nas_found <- unlist(lapply(X = data, FUN = anyNA))
  
  if(any(nas_found)) {
    stop("Variable \"", paste(valid_cols[nas_found], collapse = ", "), "\" contains missing values.")
  }
  
  # Check whether freqencies are integer ----
  not_integer <- unlist(lapply(X = data[, freq_cols], FUN = function(x) {
      any(as.integer(x)!=x)
    }
  ))
  
  if(any(not_integer)) {
    stop("Variable \"", paste(freq_cols[not_integer], collapse = ", "), "\" contains non-integer values.")
  }

  
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