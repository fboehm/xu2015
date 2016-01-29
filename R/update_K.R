
#' Calculate the allocation probability
#'
#' @param y data vector
#' @param w a scalar weight
#' @param mu a scalar class mean
#' @param kappa a scalar inverse variance
#' @export
calc_allocation_prob <- function(y, w, # w a scalar
                                 mu, # mu a scalar
                                 kappa # scalar
){
  stopifnot(length(w) == 1, length(mu) == 1, length(kappa) == 1, w <= 1, w >= 0)
  w * sqrt(kappa) * dnorm(y, mean = mu, sd = 1 / sqrt(kappa))
}


#' Calculate rho, the acceptance probability for the update K step
#'
#' @param omega_small list containing the current state of parameter vector
#' @param omega_big list containing the proposed state of parameter vector
#' @param y data vector
#' @param ind indicator denoting which cluster is affected
#' @param a hyperparameter scalar
#' @param b a scalar hyperparameter
#' @param alpha scalar parameter for dimension-matching
#' @param beta scalar parameter for dimension-matching
#' @param r scalar parameter for dimension-matching
#' @param delta scalar hyperparameter
#' @export
calc_rho <- function(y, omega_small, omega_big, ind, a, b, alpha, beta, r, delta){
  # work in one-dimension
  ##########
  # unpack omega
  K_big <- omega_big$K
  kappa_big <- omega_big$kappa[c(ind, K_big)]
  s_big <- omega_big$s # use full vector, s
  w_big <- omega_big$w[c(ind, K_big)]
  #unpack omega_prop
  kappa_small <- omega_small$kappa[ind]
  s_small <- omega_small$s
  w_small <- omega_small$w[ind]
  K_small <- omega_small$K
  # calc kappa ratio
  kappa_ratio <- (1 / gamma(a / 2)) * (kappa_big[1]) ^ (1 - a / 2) * kappa_big[2] * (b / (2 * kappa_big[2])) ^ (a / 2) * exp(- 0.5 * b * (1 / kappa_big[1] + 1 / kappa_big[2] - 1 / kappa_small)) / kappa[ind] ^ (1 - a / 2)
  # calc w ratio
  n_big <- numeric(length = 2)
  n_big[1] <- sum(s_big == ind)
  n_big[2] <- sum(s_big == K_big)
  w_ratio <- w_big[1] ^ (delta - 1 + n_big[1]) * w_big[2] ^ (delta - 1 + n_big[2]) / (w_small ^ (delta - 1 + n_new[1] + n_new[2]) * beta(delta, K_small * delta))
  # calc mu ratio
  mu_ratio <- det(calc_C(omega_big$mu)) / det(calc_C(omega_small$mu))
  # calc lik_ratio
  log_lik_ratio <- sum(dnorm(y, mean = omega_big$mu[s_big], sd = sqrt(1/omega_big$kappa)[s_big], log = TRUE) - dnorm(y, mean = omega_small$mu[s_small], sd = sqrt(1/omega_small$kappa)[s_small], log=TRUE))
  # note that we are using the full vector $y$ when calculating log lik ratio
  # alternatively, we could use only those entries of y that have s corresponding to the component that's being split.
  lik_ratio <- exp(log_lik_ratio)
  posterior_ratio <- lik_ratio * kappa_ratio * w_ratio * mu_ratio
  #### define constants
  q_K_big_d <- 0.5
  q_K_big_c <- 1 / (K_big * (K_big-1))
  qKu <- (K_small == 1) + (K_small > 1) / 2
  qKs <- 1 / K_small
  qu <- dbeta(alpha, 1, 1) * dbeta(beta, 1, 1) * dbeta(r, 2, 2)
  detJ <- (w_small ^ (3 + 1) / (w_big[1] * w_big[2]) ^ (3 / 2)) * kappa_small ^ 1.5 * (1 - r ^ 2)
  ###
  acc_ratio <- posterior_ratio * q_K_big_d * q_K_big_c * detJ / ( (K_big) * qKu * qKs * qu)
  return(acc_ratio)
}



#' Update $K$ in Gibbs sampling
#'
#' @param mu vector of class means
#' @param w vector of class weights
#' @param kappa vector of inverse variances, one entry per class
#' @param tau hyperparameter
#' @param theta hyperparameter
#' @export
update_K <- function(y, mu, w, kappa, s, tau, theta, delta){
  K <- length(mu)
  sigma <- 1 / sqrt(kappa)
  a <- 1 / (4 * tau ^ 2)
  b <- 1 / (2 * theta ^ 2)
  q_down <- 0
  if (K > 1){
    q_down <- 0.5
  }
  ### define extra parameters
  alpha <- rbeta(n = 1, 1, 1)
  beta <- rbeta(n = 1, 1, 1)
  r <- rbeta(n = 1, 2, 2)
  ### decide to split or combine
  split <- as.logical(rbinom(n = 1, size = 1, prob = 1 - q_down))
  if (split){
    omega_small <- list(K = K, mu = mu, kappa = kappa, w = w, s = s)
    ind <- sample(1:K, size = 1, replace = FALSE)
    ## edit w; make K+1 the 'new' component
    w_new <- w
    w_new[ind] <- alpha * w[ind]
    w_new[K + 1] <- (1 - alpha) * w[ind]
    ## edit mu
    mu_new <- mu
    mu_new[ind] <- mu[ind] - sqrt(w_new[K + 1] / w_new[ind]) * r / sigma[ind]
    mu_new[K + 1] <- mu[ind] + sqrt(w_new[ind] / w_new[K + 1]) * r / sigma[ind]
    ## edit kappa (sigma)
    kappa <- 1 / sigma ^ 2
    kappa_new <- kappa
    kappa_new[ind] <- beta * (1 - r) ^ 2 * (w[ind] / w_new[ind]) * kappa[ind]
    kappa_new[K + 1] <- (1 - beta) * (1 - r) ^ 2 * (w[ind] / w_new[K + 1]) * kappa[ind]
    ## Yanxun told me to see Richardson & Green 1997
    # to get the method for
    # allocation, ie, to assign values to n_new
    s_new <- s
    foo <- rbinom(n = sum(s == ind), size = 1, prob = calc_allocation_prob(y[s == ind], w_new[ind], mu_new[ind], kappa_new[ind]) / (calc_allocation_prob(y[s == ind], w_new[ind], mu_new[ind], kappa_new[ind]) + calc_allocation_prob(y[s == ind], w_new[K + 1], mu_new[K + 1], kappa_new[K + 1])))
    s_new[s == ind] <- ind * foo + (K + 1) * (1 - foo)
    ############
    omega_big <- list(K = K + 1, mu = mu_new, kappa = kappa_new, w = w_new, s = s_new)
    acc_ratio <- calc_rho(y, omega_small, omega_big, ind, a, b, alpha, beta, r)
    u <- runif(n = 1, min = 0, max = 1)
    # compare u to acceptance ratio & decide to accept or reject
    if (u < acc_ratio) {out <- list(w = w_new, mu = mu_new, kappa = kappa_new, s = s_new)} else {out <- list(w = w, mu = mu, kappa = kappa, s = s)}
  }
  ##############################################
  ### COMBINE ##################################
  else { ## combine
    indices <- sample(1:K, size=2, replace=FALSE)
    ind1 <- min(indices)
    ind2 <- max(indices)
    # edit w
    w_new <- w
    w_new[ind1]<- w[ind1] + w[ind2]
    w_new <- w_new[-ind2]
    # edit mu
    mu_new <- mu
    mu_new[ind1] <- (w[ind1] * mu[ind1] + w[ind2] * mu[ind2])/w_new[ind1]
    mu_new <- mu_new[-ind2]
    # edit kappa
    kappa_new <- kappa
    kappa_new[ind1] <- (w[ind1]/w_new[ind1])*kappa[ind1] + (w[ind2]/w_new[ind1])*kappa[ind2] + (w[ind1]*w[ind2]/(w_new[ind1])^2)*(mu[ind1] - mu[ind2])^2
    kappa_new <- kappa_new[-ind2]
    # calculate acceptance ratio
    omega_small <- list(K = K - 1, mu = mu_new, kappa = kappa_new, w = w_new, s = s_new)
    omega_big <- list(K = K, mu = mu, kappa = kappa, w = w, s = s)
    acc_ratio <- 1 / calc_rho(y, omega_small, omega_big, ind, a, b, alpha = , beta, r)
    u <- runif(n = 1, min = 0, max = 1)
    # compare u to acceptance ratio & decide to accept or reject
    if (u < acc_ratio) {out <- list(w = w_new, mu = mu_new, kappa = kappa_new, s=s_new)} else {out <- list(w = w, mu = mu, kappa = kappa, s = s)}
  }
  return(out)
}
