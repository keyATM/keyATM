#' keyATM Main function
#'
#' This is a wrapper function of \code{keyATM_fit()} and \code{keyATM_output()}.
#'
#'
#' @param keyATM_docs texts read via \code{keyATM_read()}
#' @param model keyATM model: "basic", "cov", "hmm", "lda", "ldacov" and "ldahmm"
#' @param no_keyword_topics the number of regular topics
#' @param keywords a list of keywords
#' @param model_settings a list of model specific settings
#' @param priors a list of priors of parameters
#' @param options a list of options
#' @param keep a vector of the names of elements you want to keep from \code{keyATM_fit()} output
#'
#' @return A keyATM_output object containing:
#'   \describe{
#'     \item{keyword_k}{Number of keyword topics}
#'     \item{no_keyword_topics}{Number of regular unseeded topics}
#'     \item{V}{Number of word types}
#'     \item{N}{Number of documents}
#'     \item{theta}{Normalized topic proportions for each document}
#'     \item{phi}{Normalized topic specific word generation probabilities}
#'     \item{topic_counts}{Number of tokens assigned to each topic}
#'     \item{word_counts}{Number of times each word type appears}
#'     \item{doc_lens}{Length of each document in tokens}
#'     \item{vocab}{Words in the vocabulary}
#'     \item{model_fit}{Perplexity and log-likelihood}
#'     \item{p}{Estimated p}
#'     \item{values_iter}{Organized values stored during iterations}
#'     \item{kept_values}{Output from \code{keyATM_fit()} you specified to store.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # keyATM Basic
#'   out <- keyATM(
#'                 keyATM_docs, model = "basic", no_keyword_topics = 5, keywords = keywords_list
#'                )
#'
#'   # keyATM Cov
#'   out <- keyATM(
#'                 keyATM_docs, model = "cov", no_keyword_topics = 5, keywords = keywords_list,
#'                 model_settings(covariates_data = cov)
#'                )
#'
#'   # keyATM HMM
#'   out <- keyATM(
#'                    keyATM_docs, model = "hmm", no_keyword_topics = 5, keywords = keywords_list,
#'                    model_settings(time_index = time_index_vec, num_states = 5)
#'                   )
#'
#'   # Weighted LDA
#'   out <- keyATM(
#'                 keyATM_docs, model = "lda", no_keyword_topics = 5
#'                )
#'
#'   # Weighted LDA Cov
#'   out <- keyATM(
#'                 keyATM_docs, model = "ldacov", no_keyword_topics = 5,
#'                 model_settings(covariates_data = cov)
#'                )                   
#'
#'   # Weighted LDA HMM
#'   out <- keyATM(
#'                 keyATM_docs, model = "ldahmm", no_keyword_topics = 5,
#'                 model_settings(time_index = time_index_vec, num_states = 5)
#'                )
#'
#' }
#'
#' @export
keyATM <- function(keyATM_docs, model, no_keyword_topics,
                   keywords = list(), model_settings = list(),
                   priors = list(), options = list(), keep = c())
{
  # Check type
  if (length(keep) != 0)
    check_arg_type(keep, "character")

  # Fit keyATM
  fitted <- keyATM_fit(
                       keyATM_docs, model, no_keyword_topics,
                       keywords, model_settings, priors, options
                      )

  # Get output
  out <- keyATM_output(fitted)

  # Keep some objects if specified
  if (length(keep) != 0) {
    kept_values <- list()
    use_elements <- keep[keep %in% names(fitted)]
    for(i in 1:length(use_elements)){
      kept_values[use_elements[i]]  <- fitted[use_elements[i]]
    }
    out$kept_values <- kept_values
  }

  return(out)
}


