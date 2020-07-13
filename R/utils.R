# a more than usually informative error message for handing in the
# wrong type to a function
check_arg_type <- function(arg, typename, message = NULL) {
  argname <- deparse(match.call()[['arg']])
  if (!inherits(arg, typename)) {
    if (is.null(message)) {
      stop(paste0('`', argname, '` is not a ', typename))
    } else {
      stop(message)
    }
  }
}


is.formula <- function(x) {
  inherits(x, "formula")
}


full_model_name <- function(model = c("base", "covariates", "dynamic", "label"),
                            type = c("keyATM", "lda"))
{
  model <- match.arg(model)
  type <- match.arg(type)

  if (type == "keyATM") {

    if (model == "base") {
      return("base") 
    } else if (model == "covariates") {
      return("cov") 
    } else if (model == "dynamic") {
      return("hmm") 
    } else if (model == "label") {
      return("label")
    } else {
      stop("Please select a correct model.") 
    }
  
  
  } else if (type == "lda") {

    if (model == "base") {
      return("lda") 
    } else if (model == "covariates") {
      return("ldacov") 
    } else if (model == "dynamic") {
      return("ldahmm") 
    # } else if (model == "label") {
    #   stop("Label LDA is currently not available.")
    } else {
      stop("Please select a correct model.") 
    }     
  
  } else {
    stop("Please select a correct type") 
  }

}


abb_model_name <- function(fullname)
{
  # Get abbribiation from the full name
  if (fullname %in% c("base", "lda")) {
    return("base") 
  } else if (fullname %in% c("cov", "ldacov")) {
    return("covariates") 
  } else if (fullname %in% c("hmm", "ldahmm")) {
    return("dynamic") 
  } else if (fullname %in% "label") {
    return("label")
  } else {
    stop("Invalid full model name.") 
  }

}


extract_full_model_name <- function(obj)
{
  # Get model full name from S3 class
  if ("base" %in% class(obj)) {
    return("base") 
  } else if ("cov" %in% class(obj)) {
    return("cov") 
  } else if ("hmm" %in% class(obj)) {
    return("hmm") 
  } else if ("lda" %in% class(obj)) {
    return("lda") 
  } else if ("ldacov" %in% class(obj)) {
    return("ldacov") 
  } else if ("ldahmm" %in% class(obj)) {
    return("ldahmm") 
  } else if ("label" %in% class(obj)) {
    return("label")
  }

}


rdirichlet <- function(alpha, n = 1) {
  l <- length(alpha)
  x <- matrix(stats::rgamma(l*n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x / as.vector(sm))
}


myhashmap <- function(keys, values) {
  mapped <- fastmap::fastmap(missing_default = NA)
  invisible(lapply(1:length(keys), function(x) {mapped$set(keys[x], values[x])}))
  return(mapped)
}


myhashmap_getvec <- function(mapped, keys) {
  return(vapply(keys, function(x) {mapped$get(x)}, integer(1), USE.NAMES = FALSE))
}


myhashmap_keyint <- function(keys, values) {
  mapped <- fastmap::fastmap(missing_default = NA)
  keys <- as.character(keys)  # key should be a string
  invisible(lapply(1:length(keys), function(x) {mapped$set(keys[x], values[x])}))
  return(mapped)
}


myhashmap_getvec_keyint <- function(mapped, keys) {
  keys <- as.character(keys) # key should be a string
  return(unlist(lapply(keys, function(x) {mapped$get(x)}), use.names = FALSE, recursive = FALSE))
}

standardize <- function(x) {return((x - mean(x)) / stats::sd(x))}

covariates_standardize <- function(data, type, cov_formula = NULL) {
  if (is.null(cov_formula)) {
    warning("`covariates_formula` is not provided. keyATM uses the matrix as it is.", immediate. = TRUE) 
    return(as.matrix(data))
  } else if (is.formula(cov_formula)) {
    message("Convert covariates data using `model_settings$covariates_formula`.") 
    covariates_data_use <- stats::model.matrix(cov_formula,
                                               as.data.frame(data))
  }

  if (type == "none")
    return(covariates_data_use)

  colnames_keep <- colnames(covariates_data_use)

  if (type == "all") {
    if ("(Intercept)" %in% colnames(covariates_data_use)) {
      standardize_cols <- colnames_keep[-which(colnames_keep == "(Intercept)")]
    } else {
      standardize_cols <- colnames_keep
    }
  }

  if (type == "non-factor") {
    # Ignore columns created from the factor 
    factor_cols <- names(Filter(is.factor, data))
    standardize_cols <- colnames_keep[!grepl(paste(c("^\\(Intercept\\)", paste0("^", factor_cols)), collapse = "|"), colnames_keep)]
  }

  if (length(standardize_cols) == 0) {
    return(covariates_data_use) 
  } else {
    covariates_data_use <- sapply(colnames_keep, function(col) {
                                  if (!col %in% standardize_cols)
                                    return(as.vector(covariates_data_use[, col]))
                                  return(standardize(as.vector(covariates_data_use[, col])))
                           })
    return(covariates_data_use)
  }
}

#' Convert a quanteda dictionary to keywords
#'
#' This function converts or reads a dictionary object from quanteda to a named list. "Glob"-style wildcard expressions (e.g. politic*) are resolved based on the available terms in your texts.
#' @param file file identifier for a foreign dictionary, e.g. path to a dictionary in YAML or LIWC format
#' @param docs texts read via [keyATM_read()]
#' @param dictionary a quanteda dictionary object, ignore if file is not NULL
#' @param split boolean, if multi-word terms be seperated, e.g. "air force" splits into "air" and "force".
#' @param ... additional parameters for [quanteda::dictionary()]
#' @return a named list which can be used as keywords for e.g. [keyATM()]
#' @seealso \code{\link{dictionary}}
#' @export
#' @examples
#' \dontrun{
#'   library(keyATM)
#'   library(quanteda)
#'   ## using the moral foundation dictiionary example from quanteda
#'   dictfile <- tempfile()
#'   download.file("http://bit.ly/37cV95h", dictfile)
#'   data(keyATM_data_bills)
#'   bills_dfm <- keyATM_data_bills$doc_dfm
#'   keyATM_docs <- keyATM_read(bills_dfm)
#'   read_keywords(file = dictfile, docs = keyATM_docs, format = "LIWC")
#' }
read_keywords <- function(file = NULL, docs = NULL, dictionary = NULL, split = TRUE, ...) {
    glob_s <- function(glob, wd_names) {
        stringr::str_subset(wd_names, utils::glob2rx(glob))
    }
    glob_select <- function(glob, wd_names, split, separator) {
        if (split) {
            globs <- stringr::str_split(glob, separator)[[1]]
            return(unlist(lapply(globs, glob_s, wd_names)))
        } else {
            return(glob_s(glob, wd_names))
        }
    }
    resolve_glob <- function(dict_slot, wd_names, split, separator) {
        if (is.list(dict_slot[[1]])) {
            stop("Only one-level dictionary is supported.")
        }
        unlist(lapply(dict_slot[[1]], glob_select, wd_names = wd_names, split = split, separator = separator))
    }
    if (is.null(file) & is.null(dictionary)) {
        stop("Both file and dictionary cannot be NULL. Please provide at least one of them.")
    }
    if (!is.null(file)) {
        dictionary <- quanteda::dictionary(file = file, ...)
    }
    stopifnot(class(dictionary) == "dictionary2")
    stopifnot(dictionary@meta$object$valuetype == "glob")
    if (is.null(docs)) {
        warning("'Glob'-style wildcards are not resolved.")
        return(as.list(dictionary))
    }
    wd_names <- unique(unlist(docs, use.names = FALSE, recursive = FALSE))
    res <- lapply(dictionary@.Data, resolve_glob, wd_names = wd_names, split = split, separator = dictionary@meta$object$separator)
    names(res) <- names(as.list(dictionary))
    return(res)
}
