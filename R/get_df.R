# function to generate data frame to iterate over
get_df <- function(ICC0, m0, n0, cv0, d0){

  tot <- length(ICC0)*length(m0)*length(n0)*length(d0)*length(cv0)

  ICC <- rep(ICC0, each = tot/length(ICC0))
  m <- rep(m0, each = tot/(length(ICC0)*length(m0)))
  n <- rep(n0, each = tot/(length(ICC0)*length(m0)*length(n0)))
  d <- rep(d0, each = tot/(length(ICC0)*length(m0)*length(n0)*length(d0)))
  cv <- rep(cv0, each = tot/(length(ICC0)*length(m0)*length(n0)*length(d0)*length(cv0)))

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
