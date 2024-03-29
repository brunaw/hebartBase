% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hebart.R
\name{hebart}
\alias{hebart}
\title{Hierarchical Embedded Bayesian Additive Regression Trees}
\usage{
hebart(
  formula,
  data,
  group_variable,
  num_trees = 2,
  control = list(node_min_size = 5),
  priors = list(alpha = 0.95, beta = 2, nu = 2, lambda = 0.1, tau_mu = 16 * num_trees,
    shape_sigma_phi = 0.5, scale_sigma_phi = 1, sample_sigma_phi = TRUE),
  inits = list(tau = 1, tau_phi = 1, sigma_phi = 1),
  MCMC = list(iter = 250, burn = 50, thin = 1, sigma_phi_sd = 2)
)
}
\arguments{
\item{formula}{The model formula}

\item{data}{The data to be used in the modeling}

\item{group_variable}{The name of the grouping variable}

\item{num_trees}{The number of trees (P)}

\item{control}{A list with control settings}

\item{priors}{A list with prior hyperparameters as defined by the model}

\item{inits}{A list with initial values for parameters}

\item{MCMC}{A list with MCMC parameters}
}
\value{
A list containing:
 Everything
}
\description{
This function runs a hebart model and returns the tree and
other results obtained in the last iteration of the MCMC
}
\details{
Priors used ----------------------------------------------------------
y_{ij} ~ Normal(m_j, tau^-1)
tau    ~ Gamma(nu/2, nu*lambda/2)
mu     ~ Normal(0, tau_mu^-1)
phi    ~ Normal(mu, sigma_phi^2 / T)
----------------------------------------------------------------------
}
\author{
Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
}
