#' Power calculations for simple cluster randomized trials
#'
#' Compute the power of a simple cluster randomized trial, or determine
#'   parameters to obtain a target power.
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param d The standardized effect size, or a vector of cluster sizes with
#'   length equal to M.
#' @param ICC The intra-class correlation.
#' @param M The total number of clusters. It should be even and greater than 2.
#' @param n The mean cluster size, or a vector of cluster sizes.
#' @param type The type of design effect, either 'standard' or based on
#'   the coefficient of variation 'cv'.
#' @param cv The coefficient of variation. Used when 'type' is 'cv'.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.

power.crt.test <- function(alpha = 0.05, power = 0.80,
                           M = NULL, n = NULL,
                           d = 0.20, ICC = NULL,
                           type = "standard", cv = NULL,
                           tol = .Machine$double.eps^0.25){

  # check to see that exactly one of the paramaters not specified
  num_null <- sum(sapply(list(alpha, power, M, n, d, ICC), is.null))
  if (num_null != 1) {
    stop("Exactly one of 'alpha', 'power', 'M', 'n', 'd', and 'ICC' must be NULL")
  }

  # if n is a vector of cluster sizes, calculate mean cluster size and cv
  if (length(n) > 1) {
    if (length(n) != M) {
      stop("length(n) is not equal to M. Enter a vector of the correct length,
           or enter one number for mean cluster size.")
    }
    if (type != "cv") {
      warning("length(n) > 1, so 'type' will be set to 'cv'")
      type <- "cv" # force type to be cv
    }
    n_mean <- mean(n) # find mean cluster size
    n_sd <- sd(n) # find sd of cluster sizes
    cv <- n_sd/n_mean # calculate cv
    n <- n_mean # set n to mean cluster size
  }

  # cv value needs to be present if type = cv
  if (type == "cv" & is.null(cv)) {
    stop("'cv' must be specified")
  }

  # function to get design effect
  getDEFF <- function(type, n, ICC, cv = NULL) {
    switch(type,
           "standard" = 1 + (n - 1)*ICC,
           "cv" = 1 + (((cv^2 + 1)*n) - 1)*ICC)
  }

  # create call to evaluate power
  if (is.null(n) | is.null(ICC)) {
    # if n or ICC is null, DEFF gets updated, so define DEFF inside call
    p.body <- quote({
      DEFF <- getDEFF(type, n, ICC, cv)
      qu <- qt(alpha/2, M - 2, lower.tail = FALSE)
      ncp <- sqrt((M/2)*n/(2*DEFF)) * d
      pt(qu, M - 2, ncp, lower.tail = FALSE) +
        pt(-qu, M - 2, ncp, lower.tail = TRUE)
    })
  } else {
    # DEFF doesn't depend on alpha, power, m, or d, so define outside call
    DEFF <- getDEFF(type, n, ICC, cv)
    p.body <- quote({
      qu <- qt(alpha/2, M - 2, lower.tail = FALSE)
      ncp <- sqrt((M/2)*n/(2*DEFF)) * d
      pt(qu, M - 2, ncp, lower.tail = FALSE) +
        pt(-qu, M - 2, ncp, lower.tail = TRUE)
    })
  }

  # calculate alpha
  if (is.null(alpha)) {
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     interval = c(1e-10, 1 - 1e-10),
                     tol = tol, extendInt = "yes")$root
    return(alpha)
  }

  # calculate power
  if (is.null(power)) {
    power <- eval(p.body)
    return(power)
  }

  # calculate d
  if (is.null(d)) {
    d <- uniroot(function(d) eval(p.body) - power,
                 interval = c(1e-07, 1e+07),
                 tol = tol, extendInt = "upX")$root
    return(d)
  }

  # calculate M
  if (is.null(M)) {
    M <- uniroot(function(M) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
    return(M)
  }

  # calculate n
  if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
    return(n)
  }

  # calculate ICC
  if (is.null(ICC)) {
    ICC <- uniroot(function(ICC) eval(p.body) - power,
                   interval = c(1e-7, 1e+07),
                   tol = tol, extendInt = "downX")$root
    return(ICC)
  }

}
