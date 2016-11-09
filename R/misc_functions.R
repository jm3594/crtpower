#------------------------------------------------------------------------------

# Miscellaneous functions

#------------------------------------------------------------------------------

# getClustSize

# This function generates the cluster sizes for the specified parameters:
# M is number of clusters
# n is mean cluster size
# a is min of unif distribution
# b is max of unif distribution
getClustSize <- function(M, n, a, b = n*2 - a) {
  round(runif(M, a, b))
}

#------------------------------------------------------------------------------

# completely unnecessary formatting function!!!
fix_minute <- function(x) {
  # Pass in minute(some_time)
  # if minute is less than 10, function will paste a 0 in tens place
  # looks better
  # style is everything
  ifelse(x < 10, paste(0, x, sep = ""), x)
}

#------------------------------------------------------------------------------
