#' Update w in Gibbs sampling
#' 
#' @param delta a K-length vector of hyperparameters
#' @param s a vector with one entry per observed data point of cluster assignments
#' @return w a weight vector of length $K$
#' @examples
#' update_w(c(1,1,1), c(1,1,2,2,3,3))
#' 
#' @export

update_w <- function(delta, s){
  # delta is a K-vector
  K <- length(delta)
  n <- numeric(length = K)
  # n is a length-K vector
  for (k in 1:K){
    n[k] <- sum(s == k)
  }
  out <- gtools::rdirichlet(n = 1, alpha = delta + n)
  return(out)
}
