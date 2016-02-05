#' Update tau in Gibbs sampling
#'
#' @export


update_tau <- function(theta, tau, mu, a2 = 0, b2 = 50){
  a <- 1/(4*tau^2)
  b <- 1/(2*theta^2)
  tau_prop <- theta + rnorm(n=1, mean=0, sd = 0.1)
  #tau_prop <- abs(theta + rnorm(n=1, mean=0, sd = 0.5))
  # note use of absolute value above to ensure that tau_prop is positive
  a_prop <- 1/(4*tau_prop^2)
  c <- sqrt(a^2 + 2*a*b)
  c_prop <- sqrt(a_prop^2 + 2*a_prop*b)
  # define the multiplier, mult
  mult <- b/(a + b + c)
  mult_prop <- b/(a_prop + b + c_prop)
  exponent <- 0:19
  lambda <- sqrt(2*a / (a + b + c)) * mult^exponent # vector of eigvenvalues
  lambda_prop <- sqrt(2*a_prop / (a_prop + b + c_prop))*mult_prop^exponent # vector of 20 eigenvalues
  prior_ratio <- dnorm(tau_prop, mean = a2, sd = b2) / dnorm(theta, mean = a2, sd = b2)
  e_ratio <- prod(lambda + 1) / prod(lambda_prop + 1)
  C <- calc_C(mu, theta, tau)
  C_prop <- calc_C(mu, theta, tau_prop)
  det_ratio <- det(C_prop)/det(C)
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
