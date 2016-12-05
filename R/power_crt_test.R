#' Power calculations for simple cluster randomized trials
#'
#' Compute the power of a simple cluster randomized trial, or determine
#'   parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{d}, \code{ICC}, \code{m},
#'   \code{n}, and \code{cv} must be passed as \code{NULL}. Note that
#'   \code{alpha}, \code{power}, \code{d}, and \code{cv} have non-\code{NULL}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NULL}.
#'
#' @section Note:
#'   'uniroot' is used to solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), but this function is based on
#'   work by Peter Dalgaard (power.t.test) and Stephane Champely (pwr.t.test).
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param d The standardized effect size.
#' @param ICC The intra-class correlation.
#' @param m The number of clusters per condition. It should be greater than 1.
#' @param n The mean cluster size, or a vector of cluster sizes with
#'   length equal to twice \code{m}.
#' @param cv The coefficient of variation. When \code{cv} = 0, the cluster all
#'   have the same size.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.

power.crt.test <- function(alpha = 0.05, power = 0.80,
                           m = NULL, n = NULL,
                           d = 0.20, ICC = NULL,
                           cv = 0,
                           tol = .Machine$double.eps^0.25){

  # check to see that exactly one of the paramaters not specified
  num_null <- sum(sapply(list(alpha, power, m, n, d, ICC, cv), is.null))
  if (num_null != 1) {
    stop("Exactly one of 'alpha', 'power', 'm', 'n', 'd', 'ICC', and 'cv' must be NULL")
  }

  # if n is a vector of cluster sizes, calculate mean cluster size and cv
  if (length(n) > 1) {
    if (length(n) != 2*m) {
      stop("length(n) is not equal to 2*m. Enter a vector of the correct length,
           or enter one number for mean cluster size.")
    }
    n_mean <- mean(n) # find mean cluster size
    n_sd <- sd(n) # find sd of cluster sizes
    cv <- n_sd/n_mean # calculate cv
    n <- n_mean # set n to mean cluster size
  }

  # create call to evaluate power
  if (is.null(n) | is.null(ICC)) {
    # if n or ICC is null, DEFF gets updated, so define DEFF inside call
    p.body <- quote({
      DEFF <- getDEFF(n, ICC, cv)
      qu <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)
      ncp <- sqrt(m*n/(2*DEFF)) * d
      pt(qu, 2*(m - 1), ncp, lower.tail = FALSE) +
        pt(-qu, 2*(m - 1), ncp, lower.tail = TRUE)
    })
  } else {
    # DEFF doesn't depend on alpha, power, m, or d, so define outside call
    DEFF <- getDEFF(n, ICC, cv)
    p.body <- quote({
      qu <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)
      ncp <- sqrt(m*n/(2*DEFF)) * d
      pt(qu, 2*(m - 1), ncp, lower.tail = FALSE) +
        pt(-qu, 2*(m - 1), ncp, lower.tail = TRUE)
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

  # calculate m
  if (is.null(m)) {
    m <- uniroot(function(m) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
    return(m)
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

  # calculate cv
  if (is.null(cv)) {
    cv <- uniroot(function(cv) eval(p.body) - power,
                  interval = c(1e-7, 1e+07),
                  tol = tol, extendInt = "downX")$root
    return(cv)
  }

}
