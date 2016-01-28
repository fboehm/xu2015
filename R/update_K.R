
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
  w * sqrt(kappa) * dnorm(y, mean = mu, sd = 1 / sqrt(kappa))
}


#' Calculate rho, the acceptance probability for the update K step
#'
#' @param omega_small list containing the current state of parameter vector
#' @param omega_big list containing the proposed state of parameter vector
#' @export
calc_rho <- function(y, omega_small, omega_big, a, b){
  # omega_small (AND omega_big) consists of ONLY the relevant entry (ies) of each vector:
  # mu, kappa, w;
  # and the full vector of s
  # and the integer K
  # hence, each entry in omega is a 1-vector or a 2-vector when we
  # work in one-dimension
  ##########
  # unpack omega
  kappa_big <- omega_big$kappa
  s_big <- omega_big$s
  mu_big <- omega_big$mu
  w_big <- omega_big$w
  K_big <- omega_big$K
  mu_big_full <- omega_big$mu_full
  #unpack omega_prop
  kappa_small <- omega_small$kappa
  s_small <- omega_small$s
  mu_small <- omega_small$mu
  w_small <- omega_small$w
  K_small <- omega_small$K
  mu_small_full <- omega_small$mu_full
  # calc kappa ratio
  kappa_ratio <- (1 / gamma(a / 2)) * (kappa_big[1]) ^ (1 - a / 2) * kappa_big[2] * (b / (2 * kappa_big[2])) ^ (a / 2) * exp(- 0.5 * b * (1 / kappa_big[1] + 1 / kappa_big[2] - 1 / kappa_small)) / kappa[ind] ^ (1 - a / 2)
  # calc w ratio
  w_ratio <- w_big[1] ^ (delta - 1 + n_big[1]) * w_big[2] ^ (delta - 1 + n_big[2]) / (w_small ^ (delta - 1 + n_new[1] + n_new[2]) * beta(delta, K_small * delta))
  # calc mu ratio
  mu_ratio <- det(calc_C(mu_big_full)) / det(calc_C(mu_small_full))

}



#' Update $K$ in Gibbs sampling
#'
#' @param mu vector of class means
#' @param w vector of class weights
#' @param kappa vector of inverse variances, one entry per class
#' @param tau hyperparameter
#' @param theta hyperparameter
#' @export
update_K <- function(y, mu, w, kappa, tau, theta, s){
  K <- length(mu)
  sigma <- 1 / sqrt(kappa)
  q_down <- 0
  if (K > 1){
    q_down <- 0.5
  }
  split <- as.logical(rbinom(n = 1, prob = 1 - q_down))
  if (split){
    ind <- sample(1:K, size = 1, replace = FALSE)
    a <- 1 / (4 * tau ^ 2)
    b <- 1 / (2 * theta ^ 2)
    ### define extra parameters
    alpha <- rbeta(n = 1, 1, 1)
    beta <- rbeta(n = 1, 1, 1)
    r <- rbeta(n = 1, 2, 2)
    ### define new parameters
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
    ### calculate ratios
    kappa_ratio <- (1 / gamma(a / 2)) * (kappa_new[ind]) ^ (1 - a / 2) * kappa_new[K + 1] * (b / (2 * kappa[K + 1])) ^ (a / 2) * exp(- 0.5 * b * (1 / kappa_new[ind] + 1 / kappa_new[K + 1] - 1 / kappa[ind])) / kappa[ind] ^ (1 - a / 2)
    ## Yanxun told me to see Richardson & Green 1997
    # to get the method for
    # allocation, ie, to assign values to n_new
    s_new <- s
    foo <- rbinom(n = sum(s == ind), size = 1, prob = calc_allocation_prob(y[s == ind], w_new[ind], mu_new[ind], kappa_new[ind]) / (calc_allocation_prob(y[s == ind], w_new[ind], mu_new[ind], kappa_new[ind]) + calc_allocation_prob(y[s == ind], w_new[K + 1], mu_new[K + 1], kappa_new[K + 1])))
    s_new[s == ind] <- ind * foo + (K + 1) * (1 - foo)
    n <- numeric(length=K)
    for (k in 1:K){
      n[k] <- sum(s == k)
    }
    n_new <- n
    n_new[ind] <- sum(s_new == ind)
    n_new[K + 1] <- sum(s_new == K + 1)
    ############
    w_ratio <- w_new[ind] ^ (delta - 1 + n_new[ind]) * w_new[K + 1] ^ (delta - 1 + n_new[K+1]) / (w[ind] ^ (delta - 1 + n_new[ind] + n_new[K + 1]) * beta(delta, K * delta))
    ############
    mu_ratio <- det(calc_C(mu_new)) / det(calc_C(mu))
    ############ WHAT ARE theta and tau for above calls to calc_C??
    ############
    log_lik_ratio <- sum(dnorm(y, mean = mu_new[s_new], sd = sqrt(1/kappa_new)[s_new], log = TRUE) - dnorm(y, mean = mu[s], sd = sqrt(1/kappa)[s], log=TRUE))
    lik_ratio <- exp(log_lik_ratio)
    posterior_ratio <- lik_ratio * kappa_ratio * w_ratio * mu_ratio
    #### define constants
    q_Kplus1d <- 0.5
    q_Kplus1c <- 1/(K*(K+1))
    qKu <- (K==1) + (K>1)/2
    qKs <- 1/K
    qu <- dbeta(alpha, 1, 1) * dbeta(beta, 1, 1) * dbeta(r, 2, 2)
    detJ <- (w[ind]^(3+1)/ (w_new[ind] * w_new[K+1])^(3/2)) * (kappa[ind])^1.5 * (1-r^2)
    ###
    acc_ratio <- posterior_ratio * q_Kplus1d * q_Kplus1c * detJ / ((K+1) * qKu * qKs * qu)
    u <- runif(n = 1, min = 0, max = 1)
    # compare u to acceptance ratio & decide to accept or reject
    if (u < acc_ratio) {out <- list(w_new, mu_new, kappa_new)} else {out <- list(w, mu, kappa)}
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
    w_new <- wnew[-ind2]
    # edit mu
    mu_new <- mu
    mu_new[ind1] <- (w[ind1] * mu[ind1] + w[ind2] * mu[ind2])/w_new[ind1]
    mu_new <- mu_new[-ind2]
    # edit kappa
    kappa_new <- kappa
    kappa_new[ind1] <- (w[ind1]/w_new[ind1])*kappa[ind1] + (w[ind2]/w_new[ind1])*kappa[ind2] + (w[ind1]*w[ind2]/(w_new[ind1])^2)*(mu[ind1] - mu[ind2])^2
    kappa_new <- kappa_new[-ind2]


    # calculate acceptance ratio




  }
  return(list(K=K, kappa=kappa_new, mu = mu_new, w = w_new))
}
