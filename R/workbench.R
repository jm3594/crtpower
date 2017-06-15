# delete this before packagizing

alpha <- 0.05
power <- NULL
m <- 3
n <- 100
nsd <- 0
cv <- NULL
icc <- 0.01
varw <- 0.99
varb <- NULL
d <- 0.5*sqrt(varw)

power.crt.test(alpha, power, m, n, nsd, cv, d, icc, varw, varb)
