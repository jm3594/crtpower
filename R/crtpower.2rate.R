#' Power calculations for simple cluster randomized trials, count outcome
#'
#' Compute the power of a simple cluster randomized trial with a count outcome,
#' or determine parameters to obtain a target power.
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
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param d The difference in proportions under the alternative hypothesis.
#' @param r1 The rate in the control group.
#' @param icc The intraclass correlation.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.

source("R/misc_functions.R") # remove this before packagizing

crtpower.2rate <- function(alpha = 0.05, power = 0.80,
                             m = NULL, n = NULL, cv = NULL,
                             r1 = NULL, d = NULL, icc = NULL,
                             tol = .Machine$double.eps^0.25){

  if(!is.null(m) && m <= 1) {
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

  # # checking set of p1, p2, and d
  # plist <- list(p1, p2, d)
  # pnames <- c("p1","p2","d")
  # pind <- which(unlist(lapply(plist, is.null))) # find null index
  #
  # # if one of p1, p2, and d is null, calculate it
  # if (length(pind) == 1) {
  #   assign(pnames[pind], calc_p(pind, p1, p2, d))
  # }

  r2 <- d + r1

  p.body <- quote({
    DEFF <- 1 + ((cv^2 + 1)*n - 1)*icc
    sdd <- sqrt((r1 + r2)*DEFF/(m*n))
    crit <- qpois(alpha/2, lower.tail = FALSE)
    ppois(crit - d/sdd, lower.tail = FALSE) +
      ppois(-crit - d/sdd, lower.tail = TRUE)
    # pnorm(zcrit - d/sdd, lower.tail = FALSE) +
    #   pnorm(-zcrit - d/sdd, lower.tail = TRUE)
  })


  # calculate alpha
  if (is.null(alpha)) {
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     interval = c(1e-10, 1 - 1e-10),
                     tol = tol, extendInt = "yes")$root
    #return(alpha)
  }

  # calculate power
  if (is.null(power)) {
    power <- eval(p.body)
    #return(power)
  }

  # calculate m
  if (is.null(m)) {
    m <- uniroot(function(m) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
    #return(m)
  }

  # calculate r1
  if (is.null(r1)) {
    r1 <- uniroot(function(r1) eval(p.body) - power,
                  interval = c(1e-07, 0.99999999),
                  tol = tol)$root
  }

  # calculate d
  if (is.null(d)) {
    d <- uniroot(function(d) eval(p.body) - power,
                 interval = c(1e-07, 1e+07),
                 tol = tol, extendInt = "upX")$root
    #return(d)
  }

  # calculate n
  if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
    #return(n)
  }

  # calculate cv
  if (is.null(cv)) {

    cv <- uniroot(function(cv) eval(p.body) - power,
                  interval = c(1e-7, 1e+07),
                  tol = tol, extendInt = "downX")$root
    #return(cv)
  }

  # calculate icc
  # if icc is null but varw, varb not null
  if (is.null(icc)){
    icc <- uniroot(function(icc) eval(p.body) - power,
                   interval = c(1e-07, 0.9999999),
                   tol = tol, extendInt = "downX")$root
  }

  method <- "Clustered two-sample proportion power calculation"
  note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  structure(list(m = m, n = n, cv = cv, d = d, r1 = r1, icc = icc,
                 alpha = alpha, power = power,
                 note = note, method = method),
            class = "power.htest")

}
