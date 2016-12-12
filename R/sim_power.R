# function to run the desired simulation

sim_power <- function(nsim, m, n, a, d, ICC) {
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
                       tryCatch(
                         crtPowerSim(m = m[i],
                                     n = n[i],
                                     a = a[i],
                                     B1 = d[i],
                                     B0 = 0,
                                     vare = 1 - ICC[i],
                                     varb = ICC[i]
                         ),
                         warning = function(w){
                           if(grepl(bobyqa_string, as.character(w))){
                             ww[i] <<- ww[i] + 1
                           }
                         }
                       ))

    power[i] <- mean( pvals < 0.05 )

    # do binom.test to get exact CI
    # for each row, calculate the power with crtPowerSim
    # also get upper and lower limits of exact 95% CI for power
    tmp <- binom.test(power[i]*nsim,
                      n = nsim,
                      p = 0.80)
    low[i] <- tmp$conf.int[1]
    high[i] <- tmp$conf.int[2]

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
    } # end of 'if (i == length(m))'

    #------------------------------------------------------------------------------
  } # end of 'for(i in 1:length(m))`

  return(list(low = low, power = power, high = high))
}
