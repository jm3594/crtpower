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

crtpower.t.test(alpha, power, m, n, nsd, cv, d, icc, varw, varb)

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

crtpower.2p.test(alpha, power, m, n, cv, p2, d, icc)

#-----------------------------------------------------------

# testing count outcome

m <- 28
n <- 424
cv <- 0
r1 <- 0.0148
d <- 0.0044
icc <- 0
alpha <- 0.05
power <- NULL

crtpower.2r.test(alpha, power, m, n, cv, r1, d, icc)
