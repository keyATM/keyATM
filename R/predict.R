#' Predict topic proportions for the covariate keyATM
#' 
#' @param object the keyATM_output object for the covariate model.
#' @param newdata New observations which should be predicted.
#' @param transform Transorm and standardize the `newdata` with the same formula and option as `model_settings` used in [keyATM()].
#' @param burn_in integer. Burn-in period. If not specified, it is the half of samples. Default is \code{NULL}.
#' @param parallel logical. If \code{TRUE}, parallelization for speeding up. Default is \code{TRUE}.
#' @param mc.cores integer. The number of cores to use. Default is \code{NULL}.
#' @param posterior_mean logical. If \code{TRUE}, the quantity of interest to estimate is the posterior mean. Default is \code{TRUE}.
#' @param ci value of the credible interval (between 0 and 1) to be estimated. Default is \code{0.9} (90%). 
#' @param method method for computing the credible interval. The Highest Density Interval (\code{hdi}, default) or Equal-tailed Interval (\code{eti}).
#' @param point method for computing the point estimate. \code{mean} (default) or \code{median}.
#' @param label a character. Add a `label` column to the output. The default is \code{NULL} (do not add it).
#' @param raw_values a logical. Returns raw values. The default is \code{FALSE}.
#' @param ... additional arguments not used.
#' @export
predict.keyATM_output <- function(object, newdata, transform = FALSE, burn_in = NULL, parallel = TRUE, mc.cores = NULL, 
                                  posterior_mean = TRUE, ci = 0.9, method = c("hdi", "eti"), 
                                  point = c("mean", "median"), label = NULL, raw_values = FALSE, ...)
{
  method <- match.arg(method)
  point <- match.arg(point)

  if (object$model != "covariates" | !("keyATM_output" %in% class(object)))
    stop("This is not an output of covariate model")

  if (transform) {
    if (!identical(colnames(newdata), colnames(object$kept_values$model_settings$covariates_data)))
      stop("Column names in `newdata` are different from the data provided to `keyATM()` when fitting the model.")
    newdata <- covariates_standardize(newdata, type = object$kept_values$model_settings$standardize,
                                      cov_formula = object$kept_values$model_settings$covariates_formula)
  }
  
  data_used <- covariates_get(object)
  if (dim(newdata)[1] != dim(data_used)[1] | dim(newdata)[2] != dim(data_used)[2])
    stop("Dimension of the `newdata` should match with the fitted data. Check the output of `covariates_get()`")

  if (is.null(burn_in))
    burn_in <- floor(max(object$model_fit$Iteration) / 2) 

  if (parallel) {
    if (is.null(mc.cores)) {
      num_core <- max(1, parallel::detectCores(all.tests = FALSE, logical = TRUE) - 2L)
    } else {
      num_core <- mc.cores 
    }
  } else {
    num_core <- 1L
  }

  used_iter <- object$values_iter$used_iter
  used_iter <- used_iter[used_iter > burn_in]
  use_index <- which(object$values_iter$used_iter %in% used_iter)
  tnames <- rownames(object$phi)
  Lambda_iter <- object$values_iter$Lambda_iter

  if (posterior_mean) {
    obj <- do.call(dplyr::bind_rows,
                   parallel::mclapply(1:length(use_index),
                                      function(s) {
                                        Alpha <- exp(Matrix::tcrossprod(
                                                       newdata,
                                                       Lambda_iter[[use_index[s]]]
                                                     ))
                                        rowsum <- Matrix::rowSums(Alpha)
                                        thetas <- Alpha / rowsum
                                        thetas <- as.data.frame(thetas)
                                        colnames(thetas) <- tnames
                                        thetas$Iteration <- used_iter[s]
                                        return(thetas)
                                      },
                                      mc.cores = num_core
                                     )
                  ) 
  } else{
    obj <- do.call(dplyr::bind_rows,
                   parallel::mclapply(1:length(use_index),
                                      function(s) {
                                        Alpha <- exp(Matrix::tcrossprod(
                                                       newdata,
                                                       Lambda_iter[[use_index[s]]]
                                                     ))
                                        thetas <- t(apply(Alpha, 1, rdirichlet))
                                        thetas <- Matrix::colMeans(thetas)
                                        thetas <- t(as.data.frame(thetas))
                                        thetas <- as.data.frame(thetas)
                                        colnames(thetas) <- tnames
                                        thetas$Iteration <- used_iter[s]
                                        return(thetas)
                                      },
                                      mc.cores = num_core, mc.set.seed = FALSE
                                     )
                  )
  }
  if (raw_values) {
    return(obj) 
  } else {
    return(strata_doctopic_CI(obj[, 1:(ncol(obj)-1)], ci, method, point, label))
  }
}
