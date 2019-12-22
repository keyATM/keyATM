#' Keyword Assisted Topic Models
#'
#' @section Online documentation:
#' 
#' <https://keyatm.github.io/keyATM/>
#' @useDynLib keyATM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
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

