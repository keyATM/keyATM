#' topicdict
#'
#' A package for fitting topic models informed by content
#' analysis dictionaries.
#'
#' @docType package
#' @name topicdict
#'
#' @useDynLib topicdict
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("topicdict", libpath)
}
