#' Predict topic proportions for the covariate keyATM
#'
#' @param object the keyATM_output object for the covariate model.
#' @param newdata New observations which should be predicted.
#' @param transform Transorm and standardize the `newdata` with the same formula and option as `model_settings` used in [keyATM()].
#' @param burn_in integer. Burn-in period. If not specified, it is the half of samples. Default is \code{NULL}.
#' @param parallel logical. If \code{TRUE}, parallelization for speeding up. Default is \code{TRUE}. Please `plan()` before use this function.
#' @param posterior_mean logical. If \code{TRUE}, the quantity of interest to estimate is the posterior mean. Default is \code{TRUE}.
#' @param ci value of the credible interval (between 0 and 1) to be estimated. Default is \code{0.9} (90%).
#' @param method method for computing the credible interval. The Highest Density Interval (\code{hdi}, default) or Equal-tailed Interval (\code{eti}).
#' @param point method for computing the point estimate. \code{mean} (default) or \code{median}.
#' @param label a character. Add a `label` column to the output. The default is \code{NULL} (do not add it).
#' @param raw_values a logical. Returns raw values. The default is \code{FALSE}.
#' @param ... additional arguments not used.
#' @export
predict.keyATM_output <- function(object, newdata, transform = FALSE, burn_in = NULL, parallel = TRUE,
                                  posterior_mean = TRUE, ci = 0.9, method = c("hdi", "eti"),
                                  point = c("mean", "median"), label = NULL, raw_values = FALSE, ...)
{
  method <- rlang::arg_match(method)
  point <- rlang::arg_match(point)

  if (object$model != "covariates" | !("keyATM_output" %in% class(object)))
    cli::cli_abort("This is not an output of covariate model")

  if (transform) {
    if (!identical(colnames(newdata), colnames(object$kept_values$model_settings$covariates_data)))
      cli::cli_abort("Column names in `newdata` are different from the data provided to `keyATM()` when fitting the model.")
    newdata <- covariates_standardize(newdata, type = object$kept_values$model_settings$standardize,
                                      cov_formula = object$kept_values$model_settings$covariates_formula)
  }

  data_used <- covariates_get(object)
  if (dim(newdata)[1] != dim(data_used)[1] | dim(newdata)[2] != dim(data_used)[2])
    cli::cli_abort("Dimension of the `newdata` should match with the fitted data. Check the output of `covariates_get()`")

  if (is.null(burn_in))
    burn_in <- floor(max(object$model_fit$Iteration) / 2)

  used_iter <- object$values_iter$used_iter
  used_iter <- used_iter[used_iter > burn_in]
  use_index <- which(object$values_iter$used_iter %in% used_iter)
  tnames <- rownames(object$phi)
  covariates_model <- object$kept_values$model_settings$covariates_model

  if (covariates_model == "PG") {
    Lambda_iter <- object$kept_values$model_settings$PG_params$Lambda_list
    Sigma_iter <- object$kept_values$model_settings$PG_params$Sigma_list

    if (posterior_mean) {
      # Draw from the mean
      obj <- do.call(dplyr::bind_rows,
                     future.apply::future_lapply(1:length(use_index),
                        function(s) {
                          Mu <- newdata %*% Lambda_iter[[s]]
                          D <- nrow(newdata)
                          K <- ncol(Lambda_iter[[s]]) + 1
                          Phi <- Mu
                          theta_tilda <- exp(Phi) / (1 + exp(Phi))
                          thetas <- matrix(rep(0, D*K), nrow = D, ncol = K)
                          thetas <- calc_PGtheta_R(theta_tilda, thetas, D, K)
                          thetas <- as.data.frame(thetas)
                          colnames(thetas) <- tnames
                          thetas$Iteration <- used_iter[s]
                          return(thetas)
                        },
                        future.seed = TRUE)
                     )
    } else {
      # Draw from the mean
      obj <- do.call(dplyr::bind_rows,
                     future.apply::future_lapply(1:length(use_index),
                        function(s) {
                          Mu <- newdata %*% Lambda_iter[[s]]
                          Sigma <- Sigma_iter[[s]]
                          D <- nrow(newdata)
                          K <- ncol(Lambda_iter[[s]]) + 1
                          Phi <- sapply(1:D, function(d) {
                                          Phi_d <- rmvn1(mu = Mu[d, ], Sigma = Sigma)
                                        })
                          Phi <- t(Phi)
                          theta_tilda <- exp(Phi) / (1 + exp(Phi))
                          thetas <- matrix(rep(0, D*K), nrow = D, ncol = K)
                          thetas <- calc_PGtheta_R(theta_tilda, thetas, D, K)
                          thetas <- as.data.frame(thetas)
                          colnames(thetas) <- tnames
                          thetas$Iteration <- used_iter[s]
                          return(thetas)
                        },
                        future.seed = TRUE)
                     )
    }
  } else if (covariates_model == "DirMulti") {
    Lambda_iter <- object$values_iter$Lambda_iter


    if (posterior_mean) {
      # Dir-Multi Posterior Mean
      obj <- do.call(dplyr::bind_rows,
                     future.apply::future_lapply(1:length(use_index),
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
                        future.seed = TRUE
                       )
                    )
    } else{
      # Dir-Multi Posterior Draw
      obj <- do.call(dplyr::bind_rows,
                     future.apply::future_lapply(1:length(use_index),
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
                          future.seed = TRUE
                         )
                    )
    }
  }

  if (raw_values) {
    return(obj)
  } else {
    return(strata_doctopic_CI(obj[, 1:(ncol(obj)-1)], ci = ci, method, point, label))
  }
}
