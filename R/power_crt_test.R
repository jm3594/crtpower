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
#' The user can specify either the number of clusters per condition, \code{m},
#'   or the total number of clusters, \code{M}. If \code{m} is \code{NULL}
#'   and \code{M} isn't, then \code{m} will be calculated as half of \code{M}.
#'   If both \code{m} and \code{M} are \code{NULL} and the other parameters
#'   are not, the function will return \code{m}.
#'
#' @section Note:
#'   'uniroot' is used to solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param m The number of clusters per condition. It must be greater than 1.
#' @param n The mean of the cluster sizes, or a vector of cluster sizes with
#'   length equal to twice \code{m}.
#' @param nsd The standard deviation of the cluster sizes.
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param d The difference in condition means.
#' @param icc The intraclass correlation.
#' @param varw The within-cluster variation.
#' @param varb The between-cluster variation.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.

source("R/misc_functions.R") # remove this before packagizing

power.crt.test <- function(alpha = 0.05, power = 0.80, m = NULL,
                           n = NULL, nsd = NULL, cv = NULL,
                           d = NULL, icc = NULL,
                           varw = NULL, varb = NULL,
                           tol = .Machine$double.eps^0.25){

  if(!is.null(m) & m <= 1) {
    stop("'m' must be greater than 1.")
  }

  # if n is a vector of cluster sizes, calculate mean and cv of cluster sizes
  if (length(n) > 1) {
    if (!is.null(m) && length(n) != 2*m) { # use && to evaluate !is.null(m) first
      stop("length(n) is not equal to 2*m. Enter a vector of the correct length,
           or enter one number for mean cluster size.")
    }
    nsd <- sd(n) # find sd of cluster sizes
    n <- mean(n) # find mean cluster size
    cv <- nsd/n  # find coeffient of variation
  }

  # list of all inputs
  all_params <- list(alpha, power, m, n, nsd, cv, d, icc, varw, varb)

  # get number of null values
  all_null <- sum(unlist(lapply(all_params, is.null)))

  # checking set of n, nsd, and cv
  nlist <- list(n, nsd, cv)
  nnames <- c("n","nsd","cv")
  nind <- which(unlist(lapply(nlist, is.null))) # find null index
  # check to make sure that both n and cv are not NULL
  # if only one of n and cv is specified it will be assumed that the user wants other in the pair
  if (is.null(n) & is.null(cv)){
    stop("At least one of 'n' and 'cv' must not be NULL.")
  }

  # if one of n, nsd, and cv is null, calculate it
  if (length(nind) == 1) {
    assign(nnames[nind], calc_n(nind, n, nsd, cv))
    # if num_null also is 1, then return the value just calculated
    if (all_null == 1) {
      return(get(nnames[nind]))
    }
  }

  # checking set of icc, varw, varb
  icclist <- list(icc, varw, varb)
  iccnames <- c("icc","varw","varb")
  iccind <- which(unlist(lapply(icclist, is.null))) # find null index
  # check to make sure that at least two of icc, varw, and varb is not null
  # if only one of n and cv is specified it will be assumed that the user wants other in the pair
  #   the user can calculate nsd from these two if so desired
  if (length(iccind) > 1){
    stop("At least two of 'icc', 'varw', and 'varb' must not be NULL.")
  }

  # if one of icc, varw, and varb is null, calculate it
  if (length(iccind) == 1) {
    assign(iccnames[iccind], calc_icc(iccind, icc, varw, varb))
    # if num_null also is 1, then return the value just calculated
    if (all_null == 1) {
      return(get(iccnames[iccind]))
    }
  }

  # list of needed inputs
  needed_params <- list(alpha, power, m, n, cv, d, icc, varw)
  needed_null <- sum(unlist(lapply(needed_params, is.null)))
  # check to see that exactly one needed param is null
  if(needed_null != 1){
    stop("Exactly one of 'alpha', 'power', 'm', 'n', 'cv', 'd', 'icc', and 'varw' must be NULL.")
  }

  # create call to evaluate power
  if (is.null(n) | is.null(cv) | is.null(icc)) {
    # if n, cv, or icc is null, DEFF gets updated, so define DEFF inside call
    p.body <- quote({
      DEFF <- getDEFF(n, cv, icc)
      qu <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)
      ncp <- sqrt(m*n/(2*DEFF)) * d/sqrt(varw)
      pt(qu, 2*(m - 1), ncp, lower.tail = FALSE) +
        pt(-qu, 2*(m - 1), ncp, lower.tail = TRUE)
    })
  } else {
    # DEFF doesn't depend on alpha, power, m, or d, so define outside call
    DEFF <- getDEFF(n, cv, icc)
    p.body <- quote({
      qu <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)
      ncp <- sqrt(m*n/(2*DEFF)) * d/sqrt(varw)
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

  # calculate d
  if (is.null(d)) {
    d <- uniroot(function(d) eval(p.body) - power,
                 interval = c(1e-07, 1e+07),
                 tol = tol, extendInt = "upX")$root
    return(d)
  }

  # calculate icc
  if (is.null(icc)) {
    icc <- uniroot(function(icc) eval(p.body) - power,
                   interval = c(1e-07, 1e+07),
                   tol = tol, extendInt = "downX")$root
    return(icc)
  }

  # calculate varw
  if (is.null(varw)) {
    varw <- uniroot(function(varw) eval(p.body) - power,
                    interval = c(1e-07, 1e+07),
                    tol = tol, extendInt = "upX")$root
    return(varw)
  }

}
