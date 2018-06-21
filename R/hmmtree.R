#' Fit latent-class multinomial models
#'
#' Does a lot of nice stuff
#'
#' @param dataset Character.
#' @param data A \code{data.frame} containing the data.
#' @param model A model definition, typically the path to an \code{.eqn} file.
#' @param id Character. Name of the column that contains the subject identifier.
#'   If not specified, it is assumed that each row represents observations from one participant.
#' @param condition Character. Name of the column specifying a between-subjects factor.
#'   If not specified, no between-subjects comparisons are performed.
#'
#' @keywords internal



fit_lc <- function(
  dataset
  , data
  , model
  , id = NULL
  , condition = NULL
) {
  
  HMMTreeR::lc(
    model = model
    , data = data
    , nsubj = nrow(data)
    , nclass_max = getOption("MPTmultiverse")$hmmtree$max_classes
    , nruns = getOption("MPTmultiverse")$hmmtree$n.optim
    , fi =   getOption("MPTmultiverse")$hmmtree$fisher_information
  )
}