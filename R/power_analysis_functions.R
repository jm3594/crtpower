#' Calculate the power of a cluster randomized design.

#' @param alpha The level of significance.
#' @param delta The effect size.
#' @param ICC The intra-class correlation.
#' @param M The total number of clusters.
#' @param n The mean cluster size.
#' @param type The type of design effect, either standard or based on
#'   the coefficient of variation.
#' @param cv The coefficient of variation, if type is cv.
#' @return pwr The calculated power.
get_power <- function(alpha = 0.05, delta, ICC,
                      M, n, type = "standard", cv = NULL){

  # cv value needs to be present if type = cv
  if (type == "cv" & is.null(cv)) {
    stop("The cv is NULL. Specify a cv.")
  }

  # select design effect based on type
  DEFF <- switch(type,
                 "standard" = 1 + (n - 1)*ICC,
                 "cv" = 1 + (((cv^2 + 1)*n) - 1)*ICC)

  lambda <- delta/sqrt(2*DEFF/(0.5*M*n)) # 0.5*M is number of arms
  qu <- qt(1-alpha/2, M - 2)
  power <- pt(qu, M - 2, lambda, lower.tail = FALSE) +
    pt(-qu, M - 2, lambda, lower.tail = TRUE)
  return(power)
}

#' Calculate the standardized effect size, delta.

#' @param alpha The level of significance.
#' @param power The power.
#' @param ICC The intra class correlation.
#' @param M The total number of clusters.
#' @param n The mean cluster size.
#' @param type The type of design effect, either standard or based on
#'   the coefficient of variation.
#' @param cv The coefficient of variation, if type is cv.
#' @param tol A tolerance value.
#' @return delta The effect size.

# NOTE:
# If passing vectors to the various parameters of this function,
# it must be done in a for loop.

get_delta <- function(alpha = 0.05, power = 0.80, ICC = 0.001,
                      M = 6, n = 100, type = "standard", cv = NULL,
                      tol = .Machine$double.eps^0.25){

  # cv value needs to be present if type = cv
  if (type == "cv" & is.null(cv)) {
    stop("The cv is NULL. Specify a cv.")
  }

  # select DEFF based on type
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

#' Calculate the total number of clusters.

#' @param alpha The level of significance.
#' @param power The power.
#' @param ICC The intra class correlation.
#' @param delta The effect size.
#' @param n The mean cluster size.
#' @param type The type of design effect, either standard or based on
#'   the coefficient of variation.
#' @param cv The coefficient of variation, if type is cv.
#' @param tol A tolerance value.
#' @return M The total number of clusters.

get_M <- function(alpha = 0.05, power = 0.80, ICC = 0.001,
                  delta = 0.2, n = 100, type = "standard", cv = NULL,
                  tol = .Machine$double.eps^0.25){

  # cv value needs to be present if type = cv
  if (type == "cv" & is.null(cv)) {
    stop("The cv is NULL. Specify a cv.")
  }

  # select DEFF based on type
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
  M <- uniroot(function(M) eval(p.body) - power,
                   interval = c(2 + 1e-10, 1e+07),
                   tol = tol, extendInt = "upX")$root
  return(M)
}

#' Calculate the mean cluster size.

#' @param alpha The level of significance.
#' @param power The power.
#' @param ICC The intra class correlation.
#' @param delta The effect size.
#' @param M The total number of clusters.
#' @param type The type of design effect, either standard or based on
#'   the coefficient of variation.
#' @param cv The coefficient of variation, if type is cv.
#' @param tol A tolerance value.
#' @return n The mean cluster size.

get_n <- function(alpha = 0.05, power = 0.80, ICC = 0.001,
                  delta = 0.2, M = 6, type = "standard", cv = NULL,
                  tol = .Machine$double.eps^0.25){

  # cv value needs to be present if type = cv
  if (type == "cv" & is.null(cv)) {
    stop("The cv is NULL. Specify a cv.")
  }

  # creates a function call to evaluate the power in uniroot
  p.body <- quote({
    # select DEFF based on type
    DEFF <- switch(type,
                   "standard" = 1 + (n - 1)*ICC,
                   "cv" = 1 + (((cv^2 + 1)*n) - 1)*ICC)
    qu <- qt(alpha/2, M - 2, lower.tail = FALSE)
    pt(qu, M - 2, ncp = sqrt(0.5*M*n/(2*DEFF)) * delta, lower.tail = FALSE) +
      pt(-qu, M - 2, ncp = sqrt(0.5*M*n/(2*DEFF)) * delta, lower.tail = TRUE)
  })

  # uniroot scans an interval looking for the root of listed function
  n <- uniroot(function(n) eval(p.body) - power,
               interval = c(2 + 1e-10, 1e+07),
               tol = tol, extendInt = "upX")$root
  return(n)
}
