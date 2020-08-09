# Fitting multivariate generalized linear model.  
# Authors: Hua Zhou and Yiwen Zhang 
# Edited by Shusei Eshima

#' @noRd
#' @export
DMreg <- function(Y, X, maxiters = 150, init = NULL, epsilon = 1e-06, weight = NULL,
                  display = FALSE, LRT = FALSE, 
                  parallel = FALSE, cores = NULL, cl = NULL, sys = NULL,
                  regBeta = FALSE) {

  N <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)

  if (is.null(weight)) 
    weight <- rep(1, N)
  
  if (is.null(init)) {
    init <- matrix(0.1, p, d)
    Y_rowsums <- rowSums(Y)
    for (j in 1:d) {
      fit <- my_glm.fit(x = X, y = Y[, j]/Y_rowsums, family = quasibinomial(link = "logit"))
      init[, j] <- fit$coefficients
    }
  } else {
    init <- t(init) 
  } 

  try(est <- fit_DMReg(Y = Y, X = X, weight = weight, init = init, 
                       epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                       cores = cores, cl = cl, sys = sys))

  if (is.null(est)) {
    warning("Irregularity in optimization, consider using the slice sampling.")
    return(init)
  }

  return (est)  # cov \times topic
}

fit_DMReg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
                       cores, cl, sys)
{
  d <- ncol(Y)
  p <- ncol(X)
  m <- rowSums(Y)
  N <- nrow(Y)
  beta <- init
  lliter <- rep(0, maxiters)
  lliter[1] <- ddirmnCpp(Y, init, X)
  ll2 <- lliter[1]
  niter <- 1
  
  ## ----------------------------------------##
  ## Begin the main loop
  ## ----------------------------------------##
  while ((niter <= 2 || abs(ll2 - ll1)/(abs(ll1) + 1) > epsilon) & (niter < maxiters)) {
    niter <- niter + 1
    ll1 <- lliter[niter - 1]
    res_helper <- list()
    res_helper <- objfun_helper(beta, X, Y, d, p, res_helper)
    Hessian <- -res_helper$Hessian
    dl <- -res_helper$dl
    tmpvector <- res_helper$tmpvec
    tmpmatrix <- res_helper$tmpmat
    Beta <- res_helper$Beta

    if (all(!is.na(dl)) & mean(dl^2) < 1e-04) 
      break
    temp.try <- NULL
    try(temp.try <- solve(Hessian, dl), silent = TRUE)
    if (is.null(temp.try) | any(is.nan(temp.try))) {
      ll.Newton <- NA
    } else if (is.numeric(temp.try)) {
      beta_Newton <- beta - matrix(temp.try, p, d)
      ll.Newton <- ddirmnCpp(Y, beta_Newton, X)
      ## ----------------------------------------## Half stepping
      if (is.nan(ll.Newton) || ll.Newton >= 0) {
        ll.Newton <- NA
      } else if (!is.na(ll.Newton) & ll1 >= ll.Newton) {
        for (st in 1:20) {
          beta_N <- beta - matrix(temp.try * (0.5^st), p, d)
          # llnew <- ddirmnCpp(Y, exp(X %*% beta_N))
          llnew <- ddirmnCpp(Y, beta_N, X)
          if (is.na(llnew) | is.nan(llnew) | llnew > 0) {
            next
          } else if (llnew > ll.Newton) {
            ll.Newton <- llnew
            beta_Newton <- beta_N
          }
          if (!is.na(llnew) & llnew > ll1) {
            break
          }
        }
      }
    } else {
      ll.Newton <- NA
    }
    if (is.na(ll.Newton) || ll.Newton < ll1) {
      ## ----------------------------------------## 
      ## MM update
      ## ----------------------------------------## 
      beta_MM <- beta
      weight.fit <- weight * tmpvector
      wnz <- weight.fit != 0 & weight.fit != Inf
      weight.fit <- weight.fit[wnz]
      X1 <- X[wnz, ]
      Y_new <- Beta * tmpmatrix
      Y_new <- Y_new[wnz, ]/weight.fit
      if (!parallel) {
        for (j in 1:d) {
          ## Surrogate 1 Poisson Regression
          wy <- Y_new[, j]
          wy[is.na(wy)] <- Y_new[is.na(wy), j]
          ff <- my_glm.fit(X1, Y_new[, j], weights = weight.fit, family = poisson(link = "log"), 
                        control = list(epsilon = epsilon))
          if (ff$converged) 
            beta_MM[, j] <- ff$coefficients  #par
        }
      } else {
        stop("Do not use the parallel option.") 
      }
      ll.MM <- ddirmnCpp(Y, beta_MM, X) 

      ## ----------------------------------------## 
      ## Choose the update
      ## ----------------------------------------## 
      if (is.na(ll.Newton) | (ll.MM < 0 & ll.MM > ll1)) {
        beta <- beta_MM
        lliter[niter] <- ll.MM
        ll2 <- ll.MM
      }
    } else {
      beta <- beta_Newton
      lliter[niter] <- ll.Newton
      ll2 <- ll.Newton
    }
  }
  colnames(beta) <- colnames(Y)
  return(beta) 
}


my_glm.fit <-
function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = list(), intercept = TRUE, singular.ok = TRUE) 
{  # Adapted version of stats::glm.fit()
    control <- do.call("glm.control", control)
    y <- abs(y)  # clip tiny values
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- FALSE  
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    dev.resids <- family$dev.resids
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) 
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    coefold <- NULL
    eta <- if (!is.null(etastart)) 
        etastart
    else if (!is.null(start)) 
        if (length(start) != nvars) 
            stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
              nvars, paste(deparse(xnames), collapse = ", ")), 
              domain = NA)
        else {
            coefold <- start
            offset + as.vector(if (NCOL(x) == 1L) 
              x * start
            else x %*% start)
        }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE
    for (iter in 1L:control$maxit) {
        good <- weights > 0
        mu.eta.val <- mu.eta(eta)
        good <- (weights > 0) & (mu.eta.val != 0)
        if (all(!good)) {
            conv <- FALSE
            break
        }
        z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        fit <- .Call(stats:::C_Cdqrls, x[good, , drop = FALSE] * 
            w, z * w, min(1e-07, control$epsilon/1000), check = FALSE)
        if (any(!is.finite(fit$coefficients))) {
            conv <- FALSE
            break
        }
        if (nobs < fit$rank) 
            stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
              "X matrix has rank %d, but only %d observations"), 
              fit$rank, nobs), domain = NA)
        if (!singular.ok && fit$rank < nvars) 
            stop("singular fit encountered")
        start[fit$pivot] <- fit$coefficients
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace) 
            cat("Deviance = ", dev, " Iterations - ", iter, 
              "\n", sep = "")
        boundary <- FALSE
        if (!is.finite(dev)) {
            if (is.null(coefold)) 
              stop("no valid set of coefficients has been found: please supply starting values", 
                call. = FALSE)
            warning("step size truncated due to divergence", 
              call. = FALSE)
            ii <- 1
            while (!is.finite(dev)) {
              if (ii > control$maxit) 
                stop("inner loop 1; cannot correct step size", 
                  call. = FALSE)
              ii <- ii + 1
              start <- (start + coefold)/2
              eta <- drop(x %*% start)
              mu <- linkinv(eta <- eta + offset)
              dev <- sum(dev.resids(y, mu, weights))
            }
            boundary <- TRUE
        }
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            coef <- start
            break
        }
        else {
            devold <- dev
            coef <- coefold <- start
        }
    }
  return(list(coefficients = coef, converged = conv))
}

