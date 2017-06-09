# delete this before packagizing

alpha <- 0.05
power <- NULL
m <- 6
n <- 500
nsd <- NULL
cv <- 0
icc <- 0.01
varw <- 0.99
varb <- NULL
d <- 0.2*sqrt(varw)

power.crt.test(alpha, power, m, n, nsd, cv, d, icc, varw, varb)
