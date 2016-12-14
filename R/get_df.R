#' Generate data frame to iterate over.

#' Generate a data frame containing the parameters needed to calculate
#'   power.

#' @section Note:
#' The function is currently set up for the user to find the effect size,
#'   d0 using the \code{power.crt.test} function with \code{power} = 0.80
#'   and \code{alpha} = 0.05. The user will need to modify the code in
#'   order to find other parameters.

#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com})

#' @param ICC0 A vector of intraclass correlations of interest.
#' @param m0 A vector of clusters per condition of interest.
#' @param n0 A vector of average sample sizes per cluster of interest.
#' @param cv0 A vector of coefficients of variation of interest.
#' @param d0 A vector of effect sizes of interest.
#' @return A dataframe with a number of columns equal to 5 and a number
#'  of rows equal to the product of the lengths of the input
#'  vectors.

get_df <- function(ICC0, m0, n0, cv0, d0){

  # get total number of rows for data frame
  tot <- length(ICC0)*length(m0)*length(n0)*length(d0)*length(cv0)

  ICC <- rep(ICC0, each = tot/length(ICC0))
  m <- rep(m0, each = tot/(length(ICC0)*length(m0)))
  n <- rep(n0, each = tot/(length(ICC0)*length(m0)*length(n0)))
  d <- rep(d0, each = tot/(length(ICC0)*length(m0)*length(n0)*length(d0)))
  cv <- rep(cv0, each = tot/tot)

  df <- data.frame(ICC, m, n, cv, d)

  # calculate effect size needed to get 0.80 power
  for (i in 1:nrow(df)) {
    df$d[i] <- power.crt.test(alpha = 0.05,
                              power = 0.80,
                              ICC = df$ICC[i],
                              m = df$m[i],
                              n = df$n[i],
                              cv = df$cv[i],
                              d = NULL)
  }

  # get lower limit of runif interval for average cluster size
  df$a <- df$n*(1 - sqrt(3)*df$cv)

  # create empty vectors for power (including low, high ends of 95% exact CI)
  df$low <- 0
  df$power <- 0
  df$high <- 0

  return(df)
}
