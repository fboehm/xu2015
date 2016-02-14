#' Calculate rho(omega_current, omega_proposal)
#'
#' @param omega_big parameter vector from model with bigger $K$
#' @param omega_small parameter vector from model with smaller $K$
#' @param y a data vector
#' @param theta DPP hyperparameter
#' @param tau DPP hyperparameter
#'
#' @export
calc_rho <- function(omega_small, omega_big, y, theta, tau){
  #unpack omega_small
  mu_small <- omega_small$mu
  K_small <- length(mu_small)
  w_small <- omega_small$w
  kappa_small <- omega_small$kappa
  s_small <- omega_small$s
  #unpack omega_big
  mu_big <- omega_big$mu
  K_big <- length(mu_big)
  w_big <- omega_big$w
  kappa_big <- omega_big$kappa
  s_big <- omega_big$s
  # calculate qu
  extras <- define_extra_parameters()
  qu <- dbeta(extras[1], 1, 1) * dbeta(extras[2], 1, 1) * dbeta(extras[3], 2, 2)
  # set qK_small_u
  qK_small_u <- (K_small == 1) + (K_small > 1) / 2
  # set qK_small_sj
  qK_small_s <- 1 / K_small
  # set q_K_big_d
  q_K_big_d <- 1 / 2
  # set q_K_big_c
  q_K_big_c <- 1 / (K_big * (K_big - 1))
  # ratio1: everything except 1. detJ and 2. posterior ratio
  ratio1 <- q_K_big_c * q_K_big_d / (K_big * qK_small_u * qK_small_s * qu)
  # calc det J
  detJ <- calc_detJ(w_big = w_big, w_small = w_small, kappa_big = kappa_big, kappa_small = kappa_small, r = extras[3])
  # calc w ratio
  w_ratio <- calc_w_ratio(w_big = w_big, w_small = w_small, s_big = s_big, s_small = s_small)
  # calc kappa ratio
  kappa_ratio <- calc_kappa_ratio(kappa_big = kappa_big, kappa_small = kappa_small, a = extras[1], b = extras[2])
  # calc mu ratio
  mu_ratio <- calc_mu_ratio(mu_big = mu_big, mu_small = mu_small, theta = theta, tau = tau)
  # calc likelihood ratio
  lik_ratio <- calc_lik_ratio(s_big = s_big, s_small = s_small, w_big = w_big, w_small = w_small,
                              mu_big = mu_big, mu_small = mu_small, kappa_big = kappa_big, kappa_small = kappa_small,
                              y = y)
  # posterior ratio
  out <- lik_ratio * kappa_ratio * w_ratio * mu_ratio * ratio1 * detJ
  return(out)
}

#' Calculate determinant of Jacobian
#'
#' @param w_big weight vector from model with greater $K$
#' @param w_small weight vector from model with smaller $K$
#' @param kappa_big kappa vector from model with greater $K$
#' @param kappa_small kappa vector from model with smaller $K$
#'
#' @export
calc_detJ <- function(w_big, w_small, kappa_big, kappa_small, r){
  ## first, determine w1, w1_tilde, w2_tilde
  w_int <- intersect(w_big, w_small)
  w12_tilde <- w_big[!(w_big %in% w_int)]
  w1_tilde <- w12_tilde[1]
  w2_tilde <- w12_tilde[2]
  w1 <- w_small[!(w_small %in% w_int)]
  ## second, determine kappa1 (from kappa_small)
  kappa_int <- intersect(kappa_small, kappa_big)
  kappa1 <- kappa_small[!(kappa_small %in% kappa_int)]
  # calculate determinant of J
  out <- w1 ^ (3 + 1) * kappa1 ^ (3 / 2) * (1 - r^2) / (w1_tilde * w2_tilde) ^ (3 / 2)
  return(out)
}

#' Calculate w ratio
#'
#' @param w_big w vector of weights from model with greater $K$
#' @param w_small w vector of weights from model with smaller $K$
#'
#' @export
calc_w_ratio <- function(w_big, w_small, s_big, s_small, delta = 1){
  # determine w1
  w_int <- intersect(w_big, w_small)
  w1 <- w_small[!(w_small %in% w_int)]
  # determine w1_tilde & w2_tilde
  w12_tilde <- w_big[!(w_big %in% w_int)]
  w1_tilde <- w12_tilde[1]
  w2_tilde <- w12_tilde[2]
  # determine n1_tilde and n2_tilde
  ind12 <- which(!(w_big %in% w_int))
  ind1 <- min(ind12)
  ind2 <- max(ind12)
  n1_tilde <- sum(s_big == ind1)
  n2_tilde <- sum(s_big == ind2)
  # ratio calcs
  numer <- w1_tilde ^ (delta - 1 + n1_tilde) * w2_tilde ^ (delta - 1 + n2_tilde)
  K_small <- length(w_small)
  denom <- w1 ^ (delta - 1 + n1_tilde + n2_tilde) * beta(delta, K_small * delta)
  return(numer / denom)
}

#' Calculate mu ratio
#'
#' @param mu_big mean vector for model with greater $K$
#' @param mu_small mean vector for model with smaller $K$
#' @param theta hyperparameter theta for DPP
#' @param tau hyperparameter tau for DPP
#' @export
calc_mu_ratio <- function(mu_big, mu_small, theta , tau){
  C_big <- calc_C(mu_big, theta, tau)
  C_small <- calc_C(mu_small, theta, tau)
  return(det(C_big) / det(C_small))
}


#' Calculate kappa ratio
#'
#' @param kappa_big kappa vector from the model with a greater $K$
#' @param kappa_small kappa vector from the model with smaller $K$
#' @param a hyperparameter a
#' @param b hyperparameter b
#' @export
calc_kappa_ratio <- function(kappa_big, kappa_small, a, b){
  # determine kappa1, kappa1_tilde, kappa2_tilde
  kappa_int <- intersect(kappa_big, kappa_small)
  kappa1 <- kappa_small[!(kappa_small %in% kappa_int)]
  kappa12_tilde <- kappa_big[!(kappa_big %in% kappa_int)]
  kappa1_tilde <- kappa12_tilde[1]
  kappa2_tilde <- kappa12_tilde[2]
  out <- kappa1_tilde ^ (1 - a / 2) * kappa2_tilde * (b / (kappa2_tilde * 2)) ^ (a / 2) / (kappa1 ^ (1 - a / 2) * gamma(a / 2)) * exp(-b / 2 * (1 / kappa1_tilde + 1 / kappa2_tilde - 1 / kappa1))
  return(out)
}

#' Calculate likelihood ratio
#'
#' @param s_big s for big model
#' @param s_small s for small model
#' @param w_big w for big model
#' @param w_small w for small model
#' @param mu_big mu for big model
#' @param mu_small mu for small model
#' @param kappa_big kappa for big model
#' @param kappa_small kappa for small model
#' @param y data vector
#' @export
calc_lik_ratio <- function(s_big = s_big, s_small = s_small, w_big = w_big, w_small = w_small,
                           mu_big = mu_big, mu_small = mu_small, kappa_big = kappa_big, kappa_small = kappa_small,
                           y = y){
  sd_big <- sqrt(1 / kappa_big)
  sd_small <- sqrt(1 / kappa_small)
  log_lik_big <- dnorm(y, mean = mu_big[s_big], sd = sd_big[s_big], log = TRUE)
  log_lik_small <- dnorm(y, mean = mu_small[s_small], sd = sd_small[s_small], log = TRUE)
  log_diff <- sum(log_lik_big) - sum(log_lik_small)
  return(exp(log_diff))
}
