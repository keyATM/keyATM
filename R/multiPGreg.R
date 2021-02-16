#' Run multinomial regression with Polya-Gamma augmentation
#'
#' Run multinomial regression with Polya-Gamma augmentation. There is no need to call this function directly. The keyATM Covariate internally uses this.
#'
#' @param Y Outcomes.
#' @param X Covariates.
#' @param num_topics Number of topics.
#' @param PG_params Parameters used in this function.
#' @param iter The default is \code{1}.
#' @param store_lambda The default is \code{0}.
#'
#' @export
multiPGreg <- function(Y, X, num_topics, PG_params, iter = 1, store_lambda = 0)
{  # Used in CovPG
  Phi <- PG_params$PG_Phi
  Sigma_phi <- PG_params$PG_SigmaPhi
  Lambda <- PG_params$PG_Lambda
  D <- nrow(Y)
  M <- ncol(X)
  K <- ncol(Lambda) + 1

  for (it in 1:iter) {
    # Sample Phi using Polya-Gamma
    mu_phi <- X %*% Lambda
    Phi <- multiPG_sample_Phi(Phi, Y, D, K, M, mu_phi, Sigma_phi)
  
    # Sample Lambda using Bayesian Multivariate Linear Regression
    res <- multiPG_sample_Lambda(Lambda, X, Phi, K, M) 
    Lambda <- res$Lambda
    Sigma_phi <- res$Sigma_phi
    Lambda <- Lambda
  
    if (it %% 10 == 0)
      print(it) 
  }
  
  PG_params$PG_Phi <- Phi
  PG_params$PG_SigmaPhi <- Sigma_phi
  PG_params$PG_Lambda <- Lambda
  PG_params$theta_tilda <- exp(Phi) / (1 + exp(Phi))

  if (store_lambda) {
    index <- length(PG_params$Lambda_list)
    PG_params$Lambda_list[[index + 1]] <- Lambda 
    PG_params$Sigma_list[[index + 1]] <- Sigma_phi
  }
  return(PG_params)
}


multiPG_sample_Lambda <- function(Lambda, X, Phi, K, M)
{
  Delta_0 <- diag(rep(1, M))
  tXX <- t(X) %*% X
  Delta_n <- tXX + Delta_0

  B_0 <- matrix(0, nrow = M, ncol = K-1)
  B_n <- solve(tXX + Delta_0) %*% (t(X) %*% Phi + Delta_0 %*% B_0)

  nu_0 <- K + 1
  nu_n <- nu_0 + nrow(Phi)

  V_0 <- diag(rep(1, K-1))
  XBn <- X %*% B_n
  BnB0 <- B_n - B_0
  V_n <- V_0 + t(Phi - XBn) %*% (Phi - XBn) + t(BnB0) %*% Delta_0 %*% BnB0
  
  Sigma_phi <- rinvwishart(nu = nu_n, S = V_n)

  Lambda <- matrixNormal::rmatnorm(M = B_n, U = solve(Delta_n), V = Sigma_phi)
  return(list(Lambda = Lambda, Sigma_phi = Sigma_phi))
}


multiPG_sample_Phi <- function(Phi, Y, D, K, M, mu_phi, Sigma_phi)
{
  # Cumulative sum: https://gist.github.com/rbresearch/4360311

  # Create N_dk
  apply(Y, 1,
        function(vec){
          res <- rep(NA, K)

          for (k in 1:K) {
           if (k == 1) {
             res[k] <- sum(vec) 
           } else {
             res[k] <- res[k - 1] - vec[k - 1]
             if (res[k] <= 0) {
               # can happen in simulated data 
               res[k] <- 1e-7
             }
           }
          }

          return(res)
        }) -> N_dk
  N_dk <- t(N_dk)


  # Create Kappa_dk
  Kappa_dk <- Y - N_dk/2
  Kappa_dk <- Kappa_dk[, -ncol(Kappa_dk), drop = FALSE]  # we consider K-1 categories

  # Sample Omega_dk and Phi
  Sigma_phi_inv <- solve(Sigma_phi)

  sapply(1:D,
         function(d) {

           omega <- rep(NA, K-1)
           for(k in 1:(K-1)) {
             omega[k] <- rpg(N_dk[d, k], Phi[d, k]) 
           }

           Omega_inv <- diag(1 / omega)
           Kappa_d <- Kappa_dk[d,]
           
           Sigma <- solve(diag(omega) + Sigma_phi_inv)
           mu <- Sigma %*% (Kappa_d + Sigma_phi_inv %*% mu_phi[d, ])

           Phi_d <- rmvn1(mu = mu, Sigma = Sigma)
         }) -> Phi
  Phi <- t(Phi)

  return(Phi)
}


rpg <- function(b, c)
{ # Draw from Polya-Gamma
  if (b %% 1 == 0 & (b < 20 | c == 0)) {
    return(pgdraw::pgdraw(b, c))
  } else {
    # based on Glynn et al. (2019), https://github.com/G-Lynn/DLTM/blob/master/Cpp/rpgApprox.cpp
    E_omega <- 1/(2*c) * tanh(c/2)
    V_omega <- 1.0 / (4 * c^3) * (sinh(c) - c) * (1 / cosh(c/2))^2 
    x <- stats::rnorm(n = 1, mean = b*E_omega, sd = sqrt(b*V_omega))
    return(x)
  }
}


