#' Bara et al's Abortion Debate Corpus
#'
#' A corpus of speeches from an Abortion debate in the UK House
#' of Commons, studied by Bara et al. 2007.  The unit of analysis
#' is the paragraph entry (406). Speaker and final vote metadata is included.
#'
#' @format A quanteda corpus
#' \describe{
#'   \item{speaker}{ Name of speaker}
#'   \item{vote}{ Whether they voted yes, no, or abstained}
#'   ...
#' }
#' @source Bara et al. 2007 Swiss Political Science Review
"corpus_bara_para"

#' topicdict: Topic Modeling Informed by a Content Analysis Dictionary
#'
#'
#' The \code{\link{topicdict_model}} function initializes a model from a
#' seed dictionary and
#' a file pattern specifying documents.  The \code{\link{topicdict_train}}
#' function runs a Gibbs sampler
#' on the this structure. The \code{\link{posterior}} function generates
#' useful summary statistics, e.g. theta and beta matrices.
#'
#' @docType package
#' @name topicdict
NULL

#' @useDynLib topicdict
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("topicdict", libpath)
}
