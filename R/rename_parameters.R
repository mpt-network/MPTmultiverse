#' @export
 
rename_parameters <- function(x, ...) {
  UseMethod("rename_parameters")
}

#' @export

rename_parameters.multiverseMPT <- function(x, rename = c()) {
  
  cols <- sapply(X = x, FUN = is.list)
  
  for (i in names(cols)[cols]) {
    x[[i]] <- lapply(X = x[[i]], FUN = function(x) {
      if("parameter" %in% colnames(x)) {
        x$parameter[x$parameter %in% names(rename)] <- rename[x$parameter[x$parameter %in% names(rename)]]
      }
      x
    })
  }
  x 
}