% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_K.R
\name{update_K}
\alias{update_K}
\title{Update $K$ in Gibbs sampling}
\usage{
update_K(y, mu, w, sigma, s, tau, theta, delta)
}
\arguments{
\item{y}{data vector}

\item{mu}{vector of class means}

\item{w}{vector of class weights}

\item{sigma}{vector of inverse variances, one entry per class}

\item{tau}{hyperparameter}

\item{theta}{hyperparameter}

\item{delta}{hyperparameter}
}
\description{
Update $K$ in Gibbs sampling
}

