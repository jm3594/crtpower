#' Run the desired simulation.

#' This function runs  the \code{crtPowerSim} function for each row created
#'   by the \code{get_df} function. It will also output how long it takes to
#'   do each loop. In addition, it attempts to track "bobyqa" warnings that
#'   appear as a result of using the \code{lmer} function.

#' @section NOTE:
#' The intended inputs of the function are all expected to be the same
#'   length, for example all being columns of a dataframe.

#' @param nsim The number of simulations to run for each row of the dataframe.
#' @param m The total number of clusters per condition.
#' @param n The mean cluster size.
#' @param cv The coefficient of variation. When \code{cv} = 0, the cluster all
#'   have the same size.
#' @param a The minimum of the uniform distribution used to generate cluster sizes.
#' @param d The standardized effect size.
#' @param ICC The intra-class correlation.
#' @return A list with length equal to m, n, a, d, or ICC. Each list member
#'   contains the calculated power, \code{power}, and the upper and lower bounds
#'   (\code{upper} and \code{lower}, respectively) of the exact 95% confidence
#'   interval for the power.

sim_power <- function(nsim, m, n, a, d, ICC) {

  # initialize temporary vectors for holding power and the bounds of
  # the exact 95% CI for power
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
                           if(grepl(warn_string, as.character(w))){
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
    cat("Completed loop ", i, " of ", length(m) ," in ",
        round(unclass(elapsed)[1], 2), " ",
        attr(unclass(elapsed), "units"),
        ". ", round(100*i/length(m)), "% complete. Current time is ",
        hour(then) ,":", fix_minute(minute(then)),
        ".\n", sep = "")

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
