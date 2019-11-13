# a more than usually informative error message for handing in the
# wrong type to a function
check_arg_type <- function(arg, typename, message = NULL){
  argname <- deparse(match.call()[['arg']])
  if (!inherits(arg, typename)){
    if (is.null(message)) {
      stop(paste0('`', argname, '` is not a ', typename))
    } else {
      stop(message)
    }
  }
}


is.formula <- function(x){
  inherits(x,"formula")
}


full_model_name <- function(model = c("base", "covariates", "dynamic"),
                            type = c("keyATM", "lda"))
{
  model <- measure <- match.arg(model)
  type <- measure <- match.arg(type)

  if (type == "keyATM") {

    if (model == "base") {
      return("base") 
    } else if (model == "covariates") {
      return("cov") 
    } else if (model == "dynamic") {
      return("hmm") 
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
  }else if (fullname %in% c("hmm", "ldahmm")) {
    return("dynamic") 
  }else{
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
  }

}

