% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_until_converged.R
\name{run_chains}
\alias{run_chains}
\title{Run a sampler object}
\usage{
run_chains(sampler, inits, nchains)
}
\arguments{
\item{sampler}{A named list containing \code{fun} (the function to call as the
sampler), \code{init_fun} (the function used to generate the initial
conditions), and \code{args} (a named list of arguments for the sampler
function).}

\item{inits}{Named list of initial conditions.}

\item{nchains}{Number of chains (default = 3)}
}
\description{
Run a sampler object
}
