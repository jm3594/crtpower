# delete this before packagizing

# testing continuous outcome

alpha <- 0.05
power <- 0.2233
m <- 3
n <- 100
nsd <- 0
cv <- NULL
icc <- 0.01
varw <- 0.99
varb <- NULL
d <- 0.5*sqrt(varw)

test <- crtpower.2mean(alpha, power, m, n, nsd, cv, d, icc, varw, varb)

testfun <- Vectorize(crtpower.2mean)

#----------------------------------------------------------

# testing binary outcome

m <- 20
n <- 20
cv <- 0
p1 <- 0.10
p2 <- NULL
icc <- 0.001
alpha <- 0.05
power <- 0.80
pooled <- FALSE

crtpower.2prop(alpha, power, m, n, cv, p1, p2, icc, pooled)

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
