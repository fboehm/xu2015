#' Update tau in Gibbs sampling
#'
#' @export


update_tau_gamma <- function(tau, theta, mu, shape = 3, scale = 20, exponent = 0:199){
  aa <- 1 / (4 * tau ^ 2)
  bb <- 1 / (2 * theta ^ 2)
  tau_prop <- abs(tau + rnorm(n = 1, mean = 0, sd = 1))
  #tau_prop <- abs(theta + rnorm(n=1, mean=0, sd = 0.5))
  # note use of absolute value above to ensure that tau_prop is positive
  aa_prop <- 1/(4 * tau_prop ^ 2)
  cc <- sqrt(aa ^ 2 + 2 * aa * bb)
  cc_prop <- sqrt(aa_prop ^ 2 + 2 * aa_prop * bb)
  # define the multiplier, mult
  mult <- bb / (aa + bb + cc)
  mult_prop <- bb / (aa_prop + bb + cc_prop)
  lambda <- sqrt(2 * aa / (aa + bb + cc)) * mult ^ exponent # vector of eigvenvalues
  lambda_prop <- sqrt(2*aa_prop / (aa_prop + bb + cc_prop)) * mult_prop ^ exponent # vector of 20 eigenvalues
  prior_ratio <- MCMCpack::dinvgamma(tau_prop, shape = shape, scale = scale) / MCMCpack::dinvgamma(tau, shape = shape, scale = scale)
  e_ratio <- prod(lambda + 1) / prod(lambda_prop + 1)
  C <- calc_C(mu, theta, tau)
  C_prop <- calc_C(mu, theta, tau_prop)
  det_ratio <- det(C_prop) / det(C)
  tau_acc_ratio <- prior_ratio * e_ratio * det_ratio
  u <- runif(n=1)
  if (u < tau_acc_ratio) {
    out <- tau_prop
  }
  else {
    out <- tau
  }
  return(out)
}
