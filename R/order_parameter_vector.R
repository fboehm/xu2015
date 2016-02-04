#' Order parameter vector (mu, sigma, w) according to entries in mu vector
#'
#' @param mu a mean vector
#' @param sigma a standard deviation vector
#' @param w a weight vector
#' @param s
#'
#' @export
order_parameter_vector <- function(mu, sigma, w){
  # get order
  ord <- order(mu)
  # re-assign labels to s vector
  s_new <- s
  for (i in 1:length(mu)){
    s_new[s == i] <- ord[i]
  }
  return(list(mu = mu[ord], sigma = sigma[ord], w = w[ord], s = s_new))
}
