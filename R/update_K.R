#' Calculate the allocation probability
#'
#' @param y scalar, ie, one element of data vector
#' @param w a scalar weight
#' @param mu a scalar class mean
#' @param kappa a scalar inverse variance
#' @export
calc_allocation_prob <- function(y, w, # w a scalar
                                 mu, # mu a scalar
                                 kappa # scalar
){
  log(w) +  log(sqrt(kappa)) - (y - mu) ^ 2 * kappa / 2
}




#' Update $K$ in Gibbs sampling
#'
#' @param y data vector
#' @param mu vector of class means
#' @param w vector of class weights
#' @param sigma vector of inverse variances, one entry per class
#' @param tau hyperparameter
#' @param theta hyperparameter
#' @param delta hyperparameter
#' @export
update_K <- function(y, mu, w, sigma, s, tau, theta, delta){
  kappa <- 1 / sigma^2
  a <- 1 / (4 * tau ^ 2)
  b <- 1 / (2 * theta ^ 2)
  q_combine <- 0
  if (length(mu) > 1){
    q_combine <- 0.5
  }
  # define extra parameters
  extras <- define_extra_parameters()
  ### decide to split or combine
  split <- as.logical(rbinom(n = 1, size = 1, prob = 1 - q_combine))
  ##########
  if (split){
    # define the sampling vector
    sampling_vec <- 1:(length(mu))
    ind1 <- sample(sampling_vec, size = 1, replace = FALSE)
    ind2 <- length(mu) + 1
    ## Yanxun told me to see Richardson & Green 1997
    # to get the method for
    # allocation, ie, to assign values to n_new

    # define big w
    w_big <- define_big_w(w = w, alpha = extras[1], ind1 = ind1, ind2 = ind2)
    # define big mu
    mu_big <- define_big_mu(mu = mu, w = w_big, sigma = sigma, ind1 = ind1, ind2 = ind2, r = extras[3])
    # define sigma_big... actually, kappa_big
    kappa_big <- define_big_kappa(kappa = kappa, w = w, w_new = w_big, beta = extras[2], r = extras[3], ind1 = ind1, ind2 = ind2)
    sigma_big <- 1/sqrt(kappa_big)
    # define s_big
    s_big <- define_big_s(s = s, ind1 = ind1, ind2 = ind2, y = y, w_big = w_big, mu_big = mu_big, kappa_big = kappa_big)
    ## Check that the 'new' component of mu is where it should be in the ordered mu
    if (length(mu) + 1 != length(mu_big))stop("lengths of mu and mu_big are wrong")
    mu_big[c(ind1, ind2)] -> mu_two
    mu_other <- mu_big[-c(ind1, ind2)]
    min_mu <- min(mu_two)
    max_mu <- max(mu_two)
    if (length(mu_other) == 0) good_new_mu <- TRUE
    # if mu_big has length 2, then, by definition, there is no
    # mu between the two components
    #####################################
    good_new_mu <- ! (sum(mu_other < max_mu & mu_other > min_mu) > 0)
    # in the above, good_new_mu is a logical with value TRUE if
    # the value of the new mu component is viable for acceptance
    # and FALSE otherwise
    ############
    omega_small <- list(K = length(mu), mu = mu, kappa = kappa, w = w, s = s)
    omega_big <- list(K = length(mu) + 1, mu = mu_big, kappa = kappa_big, w = w_big, s = s_big)
    #foo <- calc_rho(y, omega_small, omega_big, ind1, ind2, a, b, extras[1], extras[2], extras[3], delta = 1, theta = theta, tau = tau)
    rho <- calc_rho(omega_big, omega_small)
    u <- runif(n = 1, min = 0, max = 1)
    # 1. compare u to acceptance ratio & 2. check if good_new_mu is TRUE, then decide to accept or reject
    #print(c(good_new_mu, foo$acc_ratio))
    if (good_new_mu & u < foo$acc_ratio) {out <- list(w = w_big, mu = mu_big, kappa = kappa_big, s = s_big, ar = foo, u = u, split = split, good_new_mu = good_new_mu, indic = u < foo$acc_ratio)}
      else {out <- list(w = w, mu = mu, kappa = kappa, s = s, ar = foo, u = u, split = split, good_new_mu = good_new_mu, indic = u < foo$acc_ratio)}

  }else { ## combine
    sampling_vec <- 1:(length(mu) - 1)    # we introduce sampling_vec because there's a chance that none of the y's are assigned to some of our clusters.
    index <- sample(sampling_vec, size=1, replace=FALSE)
    ind1 <- index
    ind2 <- index + 1 # force the second index to be next to the first
    # define small w
    w_small <- define_small_w(w, ind1, ind2)
    # define small mu
    mu_small <- define_small_mu(mu = mu, w = w, w_new = w_small, ind1 = ind1, ind2 = ind2)
    # define small kappa
    kappa_small <- define_small_kappa(kappa, w, w_new = w_small, ind1, ind2, mu)
    # define small s
    s_small <- define_small_s(s, mu, ind1, ind2)
    #################
    # calculate acceptance ratio
    omega_small <- list(K = length(mu) - 1, mu = mu_small, kappa = kappa_small, w = w_small, s = s_small)
    omega_big <- list(K = length(mu), mu = mu, kappa = kappa, w = w, s = s)
    bar <- calc_rho(y, omega_small = omega_small, omega_big = omega_big, ind1, ind2, a, b, extras[1], extras[2], extras[3], delta = 1, theta = theta, tau = tau)
    acc_ratio <- 1 / bar$acc_ratio
    u <- runif(n = 1, min = 0, max = 1)
    # compare u to acceptance ratio & decide to accept or reject
    if (u < acc_ratio) {out <- list(w = w_small, mu = mu_small, kappa = kappa_small, s=s_small, ar = bar, u = u, split = split)} else {out <- list(w = w, mu = mu, kappa = kappa, s = s, ar = bar, u = u, split = split)}
  }
  return(out)
}

###########################################

#' Define three parameters, alpha, beta and r for dimension-matching purposes in RJ MCMC
#' @export
define_extra_parameters <- function(){ # for dimension-matching purposes
  ### define extra parameters
  alpha <- rbeta(n = 1, 1, 1)
  beta <- rbeta(n = 1, 1, 1)
  r <- rbeta(n = 1, 2, 2)
  return(c(alpha, beta, r))
}


###########################################
## edit w; make K+1 the 'new' component
#' Define a w weight vector that is one longer than the inputted weight vector
#' @export
define_big_w <- function(w, alpha, ind1, ind2 = length(w) + 1){
  #stopifnot(sum(w) == 1)
  #stopifnot(sum(w > 0) == length(w))
  w[ind2] <- (1 - alpha) * w[ind1]
  w[ind1] <- alpha * w[ind1]
  return(w)
}
###########################################
#' Define a mu vector
#' @export
define_big_mu <- function(mu, w, sigma, ind1, ind2 = length(w) + 1, r){
  ## edit mu
  mu[ind2] <- mu[ind1] + sqrt(w[ind1] / w[ind2]) * r / sigma[ind1]
  mu[ind1] <- mu[ind1] - sqrt(w[ind2] / w[ind1]) * r / sigma[ind1]
  return(mu)
}
###########################################

#' Define a kappa vector
#' @export
define_big_kappa <- function(kappa, w, w_new, beta, r, ind1, ind2 = length(w) + 1){
  ## edit kappa (sigma)
  kappa[ind2] <- (1 - beta) * (1 - r) ^ 2 * (w[ind1] / w_new[ind2]) * kappa[ind1]
  kappa[ind1] <- beta * (1 - r) ^ 2 * (w[ind1] / w_new[ind1]) * kappa[ind1]
  return(kappa)
}
###########################################

#' Define a s vector for the model with a longer mu vector
#' @export
define_big_s <- function(s, ind1, ind2, y, w_big, mu_big, kappa_big){
  for (i in 1:length(s)){
    if (s[i] == ind1){
      log_p1 <- calc_allocation_prob(y = y[i], w = w_big[ind1], mu = mu_big[ind1], kappa = kappa_big[ind1])
      log_p2 <- calc_allocation_prob(y = y[i], w = w_big[ind2], mu = mu_big[ind2], kappa = kappa_big[ind2])
      # calculate log denominator
      log_denom <- log(exp(log_p1) + exp(log_p2))
      log_prob <- log_p1 - log_denom
      print(c(log_denom, log_prob))
      foo_bin <- rbinom(n = 1, size = 1, prob = exp(log_prob))
      s[i] <- foo_bin * ind1 + (1 - foo_bin) * ind2
    }
  }
  return(s)
}
##########################################

#' Define a smaller (by one) weight vector
#' @export
define_small_w <- function(w, ind1, ind2){
  # edit w
  w[ind1]<- w[ind1] + w[ind2]
  w <- w[-ind2]
  return(w)
}
##########################################

#' Define a small mu vector
#' @export
define_small_mu <- function(mu, w, w_new, ind1, ind2){
  # edit mu
  mu[ind1] <- (w[ind1] * mu[ind1] + w[ind2] * mu[ind2])/w_new[ind1]
  mu <- mu[-ind2]
  return(mu)
}
##########################################

#' Define a small kappa vector
#' @export
define_small_kappa <- function(kappa, w, w_new, ind1, ind2, mu){
  kappa[ind1] <- (w[ind1] / w_new[ind1]) * kappa[ind1] + (w[ind2]/w_new[ind1])*kappa[ind2] + (w[ind1]*w[ind2]/(w_new[ind1])^2)*(mu[ind1] - mu[ind2])^2
  kappa <- kappa[-ind2]
  return(kappa)
}
##########################################

#' Define a s vector for the smaller model, ie the model with shorter mu vector
#' @export
define_small_s <- function(s, mu, ind1, ind2){
  # edit s
  s[s == ind2]<- ind1
  for (k in 1:(length(mu)-2)){
    s[s == (k+ind2)] <- k + ind2 - 1
  }
  return(s)
}

