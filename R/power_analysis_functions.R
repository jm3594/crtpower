#------------------------------------------------------------------------------

# Power Analysis Functions

# These functions are used to calculate power, effect size, sample size, etc.

#------------------------------------------------------------------------------

# getCRTPower

# This function calculates either the standard CRT Power of the CV correction

# INPUT:
# alpha is level of significance
# effect is the effect size
# M is number of clusters
# n is mean cluster size
# type is either standard DEFF or cv-based DEFF
# cv is coefficient of variation

# OUTPUT:
# pwr is the calculated power

getCRTPower <- function(alpha = 0.05, effect, ICC,
                        M, n, type = "standard", cv = NULL){

  # cv value needs to be present if type = cv
  if (type == "cv" & is.null(cv)) {
    stop("The cv is NULL. Specify a cv.")
  }

  # select design effect based on type
  DEFF <- switch(type,
                 "standard" = 1 + (n - 1)*ICC,
                 "cv" = 1 + (((cv^2 + 1)*n) - 1)*ICC)

  lambda <- effect/sqrt(2*DEFF/(0.5*M*n)) # 0.5*M is number of arms
  qu <- qt(1-alpha/2, M - 2)
  power <- pt(qu, M - 2, lambda, lower.tail = FALSE) +
    pt(-qu, M - 2, lambda, lower.tail = TRUE)
  return(power)
}

#------------------------------------------------------------------------------

# getCRTeffect

# This function calculates observable effect size

# INPUT:
# alpha is level of significance
# power is the power
# ICC is the intra class correlation
# M is number of clusters
# n is mean cluster size
# type is either standard DEFF or cv-based DEFF
# cv is coefficient of variation
# tol is a tolerance value

# OUTPUT:
# delta is the effect size

# NOTE:
# If passing vectors to the various parameters of this function,
# it must be done in a for loop.

getCRTeffect <- function(alpha = 0.05, power = 0.80, ICC = 0.001,
                         M = 6, n = 100, type = "standard", cv = NULL,
                         tol = .Machine$double.eps^0.25){

  # cv value needs to be present if type = cv
  if (type == "cv" & is.null(cv)) {
    stop("The cv is NULL. Specify a cv.")
  }

  # select design effect based on type
  DEFF <- switch(type,
                 "standard" = 1 + (n - 1)*ICC,
                 "cv" = 1 + (((cv^2 + 1)*n) - 1)*ICC)

  # creates a function call to evaluate the power in uniroot
  p.body <- quote({
    qu <- qt(alpha/2, M - 2, lower.tail = FALSE)
    pt(qu, M - 2, ncp = sqrt(0.5*M*n/(2*DEFF)) * delta, lower.tail = FALSE) +
      pt(-qu, M - 2, ncp = sqrt(0.5*M*n/(2*DEFF)) * delta, lower.tail = TRUE)
  })

  # uniroot scans an interval looking for the root of listed function
  delta <- uniroot(function(delta) eval(p.body) - power,
                   interval = c(1e-07, 1e+07),
                   tol = tol, extendInt = "upX")$root
  return(delta)
}

# change
