#------------------------------------------------------------------------------

# Miscellaneous functions

#------------------------------------------------------------------------------

# getClustSize

# This function generates the cluster sizes for the specified parameters:
# M is number of clusters
# n is mean cluster size
# a is min of unif distribution
# b is max of unif distribution
getClustSize <- function(m, n, a, b = n*2 - a) {
  round(runif(2*m, a, b))
}

#------------------------------------------------------------------------------

# minute formatting function
fix_minute <- function(x) {
  # Pass in minute(some_time)
  # if minute is less than 10, function will paste a 0 in tens place
  # looks better
  # style is everything
  ifelse(x < 10, paste(0, x, sep = ""), x)
}

#------------------------------------------------------------------------------

# function to get design effect
getDEFF <- function(n, ICC, cv) {
  1 + (((cv^2 + 1)*n) - 1)*ICC
}

#------------------------------------------------------------------------------

# function to generate data frame to iterate over
get_df <- function(ICC, m, n, cv, d){

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

#------------------------------------------------------------------------------

# function to run the desired simulation

run_sim <- function(nsim, m, n, a, d, ICC) {
  low <- numeric(length(m))
  power <- numeric(length(m))
  high <- numeric(length(m))

  for (i in 1:length(m)) {
    # record time at the beginning of the loop
    if (i == 1) {
      begin <- Sys.time()
      then <- begin # then will be updated later in the loop
      cat("Starting the first loop at ",
          hour(begin) ,":", fix_minute(minute(begin)),
          ".\n", sep = "")
    } # end of 'if (i == 1)'

    #------------------------------------------------------------------------------

    # for each ICC, M, n, and cv, initialize pvals vector
    pvals <- numeric(nsim)

    pvals <- replicate(nsim,
                       crtPowerSim(m = m[i],
                                   n = n[i],
                                   a = a[i],
                                   B1 = d[i],
                                   B0 = 0,
                                   vare = 1 - ICC[i],
                                   varb = ICC[i]))

    power[i] <- mean( pvals < 0.05 )

    #------------------------------------------------------------------------------

    # print out how long it takes each loop to run
    elapsed <- Sys.time() - then

    # update then
    then <- Sys.time()

    # unclass elapsed to get time difference
    cat("Completed loop ", i, " in ",
        round(unclass(elapsed)[1], 2), " ",
        attr(unclass(elapsed), "units"), ".\n", sep = "")

    # if this is the last loop, print ending time and duration
    if (i == length(m)) {
      total <- unclass(then - begin) # time duration
      # the last updated value for then is the ending time
      cat("Ending the last loop at ",
          hour(then), ":", fix_minute(minute(then)),
          ".\n", sep = "")
      cat("Total elapsed time: ",
          round(total[1], 2), " " ,
          attr(total,"units"), ".\n", sep = "")
    } # end of 'if (i == nrow(df))'

    #------------------------------------------------------------------------------
  } # end of 'for(i in 1:nrow(df))`

  # do binom.test to get exact CI
  # for each row, calculate the power with crtPowerSim
  # also get upper and lower limits of exact 95% CI for power
  for (i in 1:length(m)) {
    tmp <- binom.test(power[i]*nsim,
                      n = nsim,
                      p = 0.80)
    low[i] <- tmp$conf.int[1]
    high[i] <- tmp$conf.int[2]
  }

  return(list(low, power, high))
}
