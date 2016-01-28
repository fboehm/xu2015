#' Update $theta$ in Gibbs sampling


update_theta <- function(theta, tau, a1 = 0, b1 = 50){
  a <- 1/(4*tau^2)
  b <- 1/(2*theta^2)
  theta_prop <- theta + rnorm(n=1, mean=0, sd = 0.5)
  b_prop <- 1/(2*theta_prop^2)
  c <- sqrt(a^2 + 2*a*b)
  c_prop <- sqrt(a^2 + 2*a*b_prop)
  # define the multiplier, mult
  mult <- b/(a + b + c)
  mult_prop <- b_prop/(a + b_prop + c_prop)
  exponent <- 0:19
  lambda <- sqrt(2*a / (a + b + c)) * mult^exponent # vector of eigvenvalues
  lambda_prop <- sqrt(2*a / (a + b_prop + c_prop))*mult_prop^exponent # vector of 20 eigenvalues
  prior_ratio <- dnorm(theta_prop, mean = a1, sd = b1) / dnorm(theta, mean = a1, sd = b1)
  e_ratio <- prod(lambda + 1) / prod(lambda_prop + 1)
  C <- calc_C(mu, theta, tau)
  C_prop <- calc_C(mu, theta_prop, tau)
  det_ratio <- det(C_prop)/det(C)
  acc_ratio <- prior_ratio * e_ratio * det_ratio
  u <- runif(n=1)
  if (u < acc_ratio) {
    out <- theta_prop
  } 
  else {
    out <- theta
  }
  return(out)
}
