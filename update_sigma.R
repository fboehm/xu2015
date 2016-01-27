#' Update sigma in Gibbs sampling
#' 
#' @param y data vector
#' @param mu length-K vector of cluster means
#' @param s vector of cluster assignments, same length as y
#' @param a hyperparameter, a scalar
#' @param b hyperparameter, a scalar
#' @return updated sigma, a vector of cluster standard deviations
#' @examples 
#' update_sigma(rnorm(6), mu = c(0, 1, 1), s = c(1,1,2,2,3,3), a = 1, b = 1)
#' @export


update_sigma <- function(y, mu, s, a, b){
  K <- length(mu)
  n <- numeric(length=K)
  S <- numeric(length=K)
  for (k in 1:K){
    n[k] <- sum(s == k)
    S[k] <- sum((s==k)*(y - mu[k])^2)
  }
  sigmasq_inverse <- rgamma(n = K, shape = (a + n)/2, rate = (b + S)/2)
  return(sqrt(1/sigmasq_inverse))
}