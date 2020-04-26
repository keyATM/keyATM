#' keyATM with Collapsed Variational Bayes
#'
#' [Experimental feature] Fit keyATM base with Collapsed Variational Bayes
#'
#' @param docs texts read via \code{keyATM_read()}
#' @param model keyATM model: "base", "covariates", "dynamic", and "label"
#' @param no_keyword_topics the number of regular topics
#' @param keywords a list of keywords
#' @param model_settings a list of model specific settings (details are in the online documentation)
#' @param vb_options a list of settings for Variational Bayes \itemize{
#'           \item \strong{tolerance}: the default is \code{1e-4}
#'           \item \strong{initialize}: "mcmc" (default) or "random"
#'           }
#' @param priors a list of priors of parameters
#' @param options a list of options same as \code{keyATM()}. Options are used when initialization method is "mcmc".
#' @param keep a vector of the names of elements you want to keep in output
#' @return A \code{keyATM_output} object
#' @export
keyATMvb <- function(docs, model, no_keyword_topics,
                     keywords = list(), model_settings = list(), vb_options = list(),
                     priors = list(), options = list(), keep = list())
{
  message("keyATMvb is an experimental function. DO NOT USE THIS FOR THE FINAL RESULT.")
  # Check type
  if (length(keep) != 0)
    check_arg_type(keep, "character")

  model <- full_model_name(model, type="keyATM")

  # Fit keyATM
  fitted <- keyATMvb_fit(
                         docs, model, no_keyword_topics,
                         keywords, model_settings, vb_options,
                         priors, options
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

  # Add VB options
  out$vb_options <- fitted$vb_options

  return(out)
}


#' Fit a keyATM model with Collapsed Variational Bayes
#' 
#' @keywords internal
keyATMvb_fit <- function(docs, model, no_keyword_topics,
                         keywords = list(), model_settings = list(), vb_options = list(),
                         priors = list(), options = list()) 
{

}
