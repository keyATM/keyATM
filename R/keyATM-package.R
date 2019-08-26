#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL


#' keyATM: Keyword-assisted Topic Model
#'
#'
#' @docType package
#' @name keyATM
NULL


#' @useDynLib keyATM
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("keyATM", libpath)
}

