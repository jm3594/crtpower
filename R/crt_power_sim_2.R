#' Simple cluster randomized trial p-values
#'
#' Calculates the p-value under the null hypotheses of a cluster randomized
#'   trial by simulation.
#' @param M The total number of clusters. Must be even.
#' @param n The mean cluster size, or a vector of cluster sizes.
#' @param vare The variance of the measurement error.
#' @param varb The variance of the random effects.
#' @param B0 The overall intercept.
#' @param B1 The treatment effect.
#' @return The p-value under the null hypothesis.
crtPowerSim <- function(M = 6, n = 100,
                        vare = 0.95, varb = 0.05,
                        B0 = 0, B1 = 0.2) {
  sde <- sqrt(vare)
  sdb <- sqrt(varb)
  # create n_vec to hold the sample size in each cluster
  # if n is a scalar, make a copy of n for each cluster, M total
  if (length(n) == 1) {
    n_vec <- rep(n, times = M)
  } else {
    # if n is alread a vector of the correct size, use it for n_vec
    if (length(n) != M) {
      stop("length(n) is not equal to M")
    }
    n_vec <- n
  }
  N <- sum(n_vec) # total sample size
  b <- rep(rnorm(M, mean = 0, sd = sdb), times = n_vec) # random effects
  e <- rnorm(N, mean = 0, sd = sde) # random errors
  i <- rep(1:M, times = n_vec) # cluster ids
  groups <- rep(c(0,1), each = M/2) # set up groups, 0 = ctrl, 1 = trtmnt
  x <- rep(groups, times = n_vec) # assign groups
  y <- B0 + x*B1 + b + e # create response vector
  crt.mod <- lmer(y ~ x + (1|i)) # create lmm for crt

  tscore <- coef(summary(crt.mod))[2,"t value"] # extract needed tscore
  pval <- 2*(1 - pt(abs(tscore), M - 2)) # get p-value
  # for low M, p-values don't match up super well with Donner-Klar values

  return(pval)

}

