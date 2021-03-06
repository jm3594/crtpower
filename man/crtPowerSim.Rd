% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crt_power_sim.R, R/crt_power_sim_2.R
\name{crtPowerSim}
\alias{crtPowerSim}
\title{Simple cluster randomized trial p-values}
\usage{
crtPowerSim(M = 6, n = 100, vare = 0.95, varb = 0.05, B0 = 0,
  B1 = 0.2)

crtPowerSim(M = 6, n = 100, vare = 0.95, varb = 0.05, B0 = 0,
  B1 = 0.2)
}
\arguments{
\item{M}{The total number of clusters. Must be even.}

\item{n}{The mean cluster size.}

\item{vare}{The variance of the measurement error.}

\item{varb}{The variance of the random effects.}

\item{B0}{The overall intercept.}

\item{B1}{The standardized treatment effect.}

\item{a}{The minimum of the uniform distribution used to generate cluster sizes.}

\item{M}{The total number of clusters. Must be even.}

\item{n}{The mean cluster size, or a vector of cluster sizes.}

\item{vare}{The variance of the measurement error.}

\item{varb}{The variance of the random effects.}

\item{B0}{The overall intercept.}

\item{B1}{The treatment effect.}
}
\value{
The p-value under the null hypothesis.

The p-value under the null hypothesis.
}
\description{
Calculates the p-value of a cluster randomized trial by simulation.

Calculates the p-value under the null hypotheses of a cluster randomized
  trial by simulation.
}

