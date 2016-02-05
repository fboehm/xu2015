#' little c, when we're using squared exponential kernel
#'
#' @export

calc_little_c <- function(mu1, mu2, theta){#mu1 a scalar and mu2 a scalar
  stopifnot(length(theta) == 1, length(mu1) == 1, length(mu2) == 1)
  return(exp(- (mu1 - mu2)^2 / theta^2))
}

#' little q, for univariate y
#'
#' @export
calc_little_q <- function(mu, tau){
  #stopifnot(tau > 0, length(tau) == 1)
  return(dnorm(mu, mean = 0, sd = tau))
}



#' Calculate the C matrix for DPP
#'
#' @param mu vector of means
#' @param theta hyperparameter, a scalar
#' @param tau hyperparameter, a scalar
#' @export

calc_C <- function(mu, theta, tau){# mu is a numeric vector; theta & tau are scalars
  stopifnot(length(theta) == 1, length(tau) == 1)
  K <- length(mu)
  C <- matrix(NA, nrow = K, ncol = K)
  for (i in 1:K){
    for (j in 1:K){
      C[i,j]<- calc_little_q(mu[i], tau) * calc_little_c(mu[i], mu[j], theta) * calc_little_q(mu[j], tau)
    }
  }
  return(C)
}

#' Update $\mu$ in Gibbs sampling
#'
#' @param y data vector
#' @param s vector of class assignments
#' @param mu vector of class means
#' @param sigma vector of class standard deviations
#' @param tau hyperparameter
#' @export
update_mu <- function(y, s, mu, sigma, tau, theta){
  stopifnot(length(tau) == 1, length(theta) == 1, length(s) == length(s), length(mu) == length(sigma))
  K <- length(mu)
  mu_prop <- mu
  for (k in 1:K){
    eps <- rnorm(n = 1, mean = 0, sd = 0.3) # reasonable sd?
    bar <- dnorm(y[s==k], mean = mu[k], sd = sigma[k], log = TRUE)
    foo <- sum(bar)
    bar_prop <- dnorm(y[s==k], mean = mu[k] + eps, sd = sigma[k], log = TRUE)
    foo_prop <- sum(bar_prop)
    ## define C & Cprop
    C <- calc_C(mu, theta, tau)
    mu_prop[k] <- mu[k] + eps
    C_prop <- calc_C(mu_prop, theta, tau)
    b <- C[k, -k]
    b_prop <- C_prop[k, -k]
    # caution with a single mu
    if (length(mu) > 1) {bar <- C[k, k] - b %*% solve(C[-k, -k]) %*% b} else
      bar <- as.numeric(C[k,k])
    if (length(mu_prop) > 1) {bar_prop <- C_prop[k, k] - b_prop %*% solve(C_prop[-k, -k]) %*% b_prop} else
      bar_prop <- as.numeric(C_prop[k,k])
    p_ratio <- exp(foo_prop - foo) * bar_prop / bar
    u <- runif(n = 1, min = 0, max = 1)
    if ( u < p_ratio) {
      mu[k]<- mu[k] + eps
    }
  }
  return(mu)
}

