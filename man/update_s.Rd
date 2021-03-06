% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_s.R
\name{update_s}
\alias{update_s}
\title{Update s parameter in Gibbs sampling}
\usage{
update_s(y, mu, sigma, w)
}
\arguments{
\item{y}{data vector}

\item{mu}{vector of K cluster means}

\item{sigma}{vector K cluster standard deviations}

\item{w}{vector of K weights}
}
\value{
s vector of same length as y
}
\description{
Update s parameter in Gibbs sampling
}
\examples{
update_s(rnorm(n=6), mu = c(1:3), sigma = c(1,1,1), w=1:3/6)
}

