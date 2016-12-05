# Run a simulation

library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

source("misc_functions.R")
source("crt_power_sim.R")
source("power_crt_test.R")

ICC <- rep(c(0.01, 0.05), each = 108)
m <- rep(c(5,25,50), times = 6, each = 12)
n <- rep(c(10,30,100), times = 2, each = 36)
cv <- rep(seq(0,0.55, by = 0.05), times = 18)
d <- 0

df <- get_df(ICC, m, n, cv, d)

df <- head(df)

ns <- vector("list", nrow(df))

run_sim(nsim = 100, m = df$m, n = df$n, a = df$a, d = df$d, ICC = df$ICC)
