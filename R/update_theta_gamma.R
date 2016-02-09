#' Update $theta$ in Gibbs sampling
#'
#' @param theta scalar
#' @param tau scalar
#' @param mu vector of means
#' @param shape hyperparameter, prior shape for theta
#' @param scale hyperparameter, prior scale for theta
#' @export

update_theta_gamma <- function(theta, tau, mu, shape = 3, scale = 20, exponent = 0:199){
  # convert from theta & tau to a & b
  aa <- 1 / (4 * tau ^ 2)
  bb <- 1 / (2 * theta ^ 2)
  ## make a theta_prop. What is a good sd to have here?
  ## I have no intuition about what's a reasonable value for theta & tau!
  theta_prop <- abs(theta + rnorm(n = 1, mean = 0, sd = 10))
  bb_prop <- 1 / (2 * theta_prop ^ 2)
  cc <- sqrt(aa ^ 2 + 2 * aa * bb)
  cc_prop <- sqrt(aa ^ 2 + 2 * aa * bb_prop)
  # define the multiplier, mult
  mult <- bb / (aa + bb + cc)
  mult_prop <- bb_prop / (aa + bb_prop + cc_prop)
  lambda <- sqrt(2 * aa / (aa + bb + cc)) * mult ^ exponent # vector of eigvenvalues
  lambda_prop <- sqrt(2 * aa / (aa + bb_prop + cc_prop)) * mult_prop ^ exponent # vector of 20 eigenvalues
  prior_ratio <- MCMCpack::dinvgamma(theta_prop ^ 2, shape = shape, scale = scale) / MCMCpack::dinvgamma(theta ^ 2, shape = shape, scale = scale)
  e_ratio <- prod(lambda + 1) / prod(lambda_prop + 1)
  # note that lambda_prop is in denominator of e_ratio
  C <- calc_C(mu, theta, tau)
  C_prop <- calc_C(mu, theta_prop, tau)
  det_ratio <- det(C_prop) / det(C)
  acc_ratio <- prior_ratio * e_ratio * det_ratio
  u <- runif(n = 1)
  if (u < acc_ratio) {
    out <- theta_prop
  }
  else {
    out <- theta
  }
  return(out)
}
