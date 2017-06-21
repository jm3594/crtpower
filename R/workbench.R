# delete this before packagizing

# testing continuous outcome

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

#----------------------------------------------------------

# testing binary outcome

m <- 25
n <- NULL
cv <- 0
p1 <- 0.51
p2 <- 0.44
d <- 0.07
icc <- 0.02
alpha <- 0.05
power <- 0.90

power.crt2p.test(alpha, power, m, n, cv, p2, d, icc)

#-----------------------------------------------------------

# testing count outcome

m <- 40
n <- 20
cv <- 0
r1 <- 0.5
d <- 0.1
icc <- 0.02
alpha <- 0.05
power <- 0.90

power.crt2rate.test(alpha, power, m, n, cv, r1, d, icc)
