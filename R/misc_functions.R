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
