#' keyATM Main function
#'
#' Run keyATM models.
#'
#'
#' @param keyATM_docs texts read via \code{keyATM_read()}
#' @param model keyATM model: "base", "covariates", and "dynamic"
#' @param no_keyword_topics the number of regular topics
#' @param keywords a list of keywords
#' @param model_settings a list of model specific settings
#' @param priors a list of priors of parameters
#' @param options a list of options \describe{
#'      \item{seed}{a numeric value for random seed. If it is not provided, the package randomly selects a seed.}
#' }
#' @param keep a vector of the names of elements you want to keep in output
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
#'     \item{kept_values}{Outputs you specified to store.}
#'   }
#'
#' @seealso \url{https://keyatm.github.io/keyATM/articles/pkgdown_files/Options.html}
#'
#' @examples
#' \dontrun{
#'   # keyATM Base
#'   out <- keyATM(
#'                 docs, model = "base", no_keyword_topics = 5, keywords = keywords_list
#'                )
#'
#'   # keyATM Covariates
#'   out <- keyATM(
#'                 docs, model = "covariates", no_keyword_topics = 5, keywords = keywords_list,
#'                 model_settings(covariates_data = cov)
#'                )
#'
#'   # keyATM Dynamic
#'   out <- keyATM(
#'                 docs, model = "dynamic", no_keyword_topics = 5, keywords = keywords_list,
#'                 model_settings(time_index = time_index_vec, num_states = 5)
#'                )
#'
#' }
#'
#' @export
keyATM <- function(docs, model, no_keyword_topics,
                   keywords = list(), model_settings = list(),
                   priors = list(), options = list(), keep = c())
{
  # Check type
  if (length(keep) != 0)
    check_arg_type(keep, "character")

  model <- full_model_name(model, type="keyATM")

  # Fit keyATM
  fitted <- keyATM_fit(
                       docs, model, no_keyword_topics,
                       keywords, model_settings, priors, options
                      )

  # Get output
  out <- keyATM_output(fitted)

  # Keep some objects if specified
  keep <- check_arg_keep(keep, model)
  if (length(keep) != 0) {
    kept_values <- list()
    use_elements <- keep[keep %in% names(fitted)]
    for (i in 1:length(use_elements)) {
      kept_values[use_elements[i]]  <- fitted[use_elements[i]]
    }
    out$kept_values <- kept_values
  }

  # A bit of clean up
  if (fitted$options$store_theta && "stored_values" %in% keep) {
    # The same information
    out$kept_values$stored_values$Z_tables <- NULL
  }


  return(out)
}


check_arg_keep <- function(obj, model)
{
  if (model %in% c("cov", "ldacov")) {
    if (!"stored_values" %in% obj)
      obj <- c("stored_values", obj) 

    if (!"model_settings" %in% obj)
      obj <- c("model_settings", obj) 
  }

  return(obj)
}



#' Weighted LDA Main function
#'
#' Run weighted LDA models.
#'
#'
#' @param docs texts read via \code{keyATM_read()}
#' @param model Weighted LDA model: "base", "covariates", and "dynamic"
#' @param number_of_topics the number of regular topics
#' @param model_settings a list of model specific settings
#' @param priors a list of priors of parameters
#' @param options a list of options
#' @param keep a vector of the names of elements you want to keep in output
#'
#' @return A keyATM_output object containing:
#'   \describe{
#'     \item{number_of_topics}{Number of regular unseeded topics}
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
#'     \item{kept_values}{Outputs you specified to store.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Weighted LDA
#'   out <- weightedLDA(
#'                      keyATM_docs, model = "base", number_of_topics = 5
#'                     )
#'
#'   # Weighted LDA Covariates
#'   out <- weightedLDA(
#'                      keyATM_docs, model = "covariates", number_of_topics = 5,
#'                      model_settings(covariates_data = cov_matrix)
#'                     )                   
#'
#'   # Weighted LDA Dynamic
#'   out <- weightedLDA(
#'                      keyATM_docs, model = "dynamic", number_of_topics = 5,
#'                      model_settings(time_index = time_index_vec, num_states = 5)
#'                     )
#'
#' }
#'
#' @export
weightedLDA <- function(docs, model, number_of_topics,
                        model_settings = list(),
                        priors = list(), options = list(), keep = c())
{
  # Check type
  if (length(keep) != 0)
    check_arg_type(keep, "character")

  model <- full_model_name(model, type="lda")

  # Fit keyATM
  fitted <- keyATM_fit(
                       docs, model, number_of_topics,
                       keywords = list(),
                       model_settings = model_settings,
                       priors = priors,
                       options = options
                      )

  # Get output
  out <- keyATM_output(fitted)
  out$number_of_topics <- number_of_topics
  out$no_keyword_topics <- NULL
  out$keyword_k <- NULL

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


