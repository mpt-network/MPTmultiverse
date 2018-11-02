#' @keywords internal

get_eqn_categories <- function (model.filename)
{
  parse.eqn <- function(x) {
    branches <- unique(x[, 2])
    l.tree <- length(branches)
    tree <- vector("expression", l.tree)
    for (branch in 1:l.tree) {
      tree[branch] <- parse(text = paste(x[x$V2 == branches[branch],
                                           "V3"], collapse = " + "))
    }
    tree
  }
  #browser()
  tmp.in <- utils::read.table(model.filename, skip = 1, stringsAsFactors = FALSE)
  tmp.ordered <- tmp.in[order(tmp.in$V1), ]
  tmp.spl <- split(tmp.ordered, factor(tmp.ordered$V1))
  tmp.spl <- lapply(tmp.spl, function(d.f) d.f[order(d.f[, 2]), ])
  unlist(lapply(tmp.spl, function(x) unique(x$V2)))
  # model <- lapply(tmp.spl, parse.eqn)
  # names(model) <- NULL
  # model
}
