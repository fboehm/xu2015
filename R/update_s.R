#' Update s parameter in Gibbs sampling
#'
#' @param y data vector
#' @param mu vector of K cluster means
#' @param sigma vector K cluster standard deviations
#' @param w vector of K weights
#' @return s vector of same length as y
#' @examples
#' update_s(rnorm(n=6), mu = c(1:3), sigma = c(1,1,1), w=1:3/6)
#' @export

update_s <- function(y, mu, sigma, w){
  K <- length(mu)
  n_subject <- length(y)
  s <- numeric(length = n_subject)
  for (i in 1:n_subject){
    normal_lik <- dnorm(x = y[i], mean=mu, sd = sigma)
    # mu & sigma are length-K vectors, so normal_lik is a K-vector
    s[i]<- sample(x = 1:K, size = 1, prob = w * normal_lik)
    # normal_lik is a K-vector; w is a length-K vector
  }
  return(s)
}
