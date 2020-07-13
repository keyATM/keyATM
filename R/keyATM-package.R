#' Keyword Assisted Topic Models
#'
#' @description The implementation of keyATM models.
#' @importFrom Rcpp sourceCpp
#' @useDynLib keyATM, .registration = TRUE
#' @aliases NULL keyATM-package
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("keyATM", libpath)
}

