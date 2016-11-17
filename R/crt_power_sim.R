#' Simple cluster randomized trial p-values
#'
#' Calculates the p-value of a cluster randomized trial by simulation.
#' @param m The number of clusters per arm.
#' @param n The mean cluster size.
#' @param a The minimum of the uniform distribution used to generate cluster sizes.
#' @param vare The variance of the measurement error.
#' @param varb The variance of the random effects.
#' @param B0 The overall intercept.
#' @param B1 The standardized treatment effect.
#' @return The p-value under the null hypothesis.
crtPowerSim <- function(M = 10, n = 10, a = 10,
                        vare = 0.95, varb = 0.05,
                        B0 = 0, B1 = 0.2) {
  sde <- sqrt(vare)
  sdb <- sqrt(varb)
  n_vec <- round(runif(M, a, 2*n - a))
  # n_vec <- getClustSize(M, n, a) # M cluster sizes from unif w/ mean n
  N <- sum(n_vec) # total sample size
  b <- rep(rnorm(M, sd = sdb), times = n_vec) # generate random effects
  e <- rnorm(N, sd = sde) # generate random errors
  i <- rep(1:M, times = n_vec) # generate cluster ids
  groups <- rep(c(0,1), each = M/2) # set up groups, 0 = ctrl, 1 = trtmnt
  # randassign <- sample(groups) # get a permutation of groups
  x <- rep(groups, times = n_vec) # random assign groups
  y <- B0 + x*B1 + b + e # create response vector
  crt.mod <- lmer(y ~ x + (1|i)) # create lmm for crt

  tscore <- coef(summary(crt.mod))[2,"t value"] # extract needed tscore
  pval <- 2*(1 - pt(abs(tscore), M - 2)) # get p-value
  # for low M, p-values don't match up super well with Donner-Klar values

  return(pval)

}

