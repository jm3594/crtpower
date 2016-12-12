# Run a simulation

library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

setwd("~/crtpower/R")

source("misc_functions.R")
source("crt_power_sim.R")
source("power_crt_test.R")
source("get_df.R")
source("sim_power.R")

# initialize vectors to pass to get_df() function
ICC0 <- c(0.01, 0.05)
m0 <- c(5, 25, 50)
n0 <- c(10, 30, 100)
d0 <- 0
cv0 <- seq(0,0.55, by = 0.05)

# alternate set of vectors to pass to get_df() function
# ICC0 <- c(0.01, 0.05, 0.01, 0.05)
# m0 <- c(2, 5, 10, 25)
# n0 <- c(10, 20, 50, 100)
# d0 <- 0
# cv0 <- 0

# create data frame to iterate over
df <- get_df(ICC0, m0, n0, cv0, d0)

# bobyqa warning counter for each row of data frame
ww <- numeric(nrow(df))

# string to capture with using tryCatch
bobyqa_string <- "bobyqa"

nsim <- 100

# simulate power
power_list <- sim_power(nsim,
                        m = df$m, n = df$n, a = df$a,
                        d = df$d, ICC = df$ICC)

# store the output of the simulation in df
df$low <- power_list$low
df$power <- power_list$power
df$high <- power_list$high

# labels for m and n
mlabs <- c(`5` = "M = 10", `25` = "M = 50", `50` = "M = 100")
nlabs <- c(`10` = "n = 10", `30` = "n = 30", `100` = "n = 100")

# subset by ICC
df1 <- filter(df, ICC == 0.01)
df2 <- filter(df, ICC == 0.05)

# generate plots
(dfplot1 <- ggplot(df1, aes(x = cv, y = power)) + ylim(0,1) +
    geom_abline(aes(slope = 0, intercept = 0.8, color = "red")) +
    geom_linerange(aes(ymin = low, ymax = high), size = 3) +
    facet_grid(n ~ m, labeller = labeller(m = as_labeller(mlabs),
                                              n = as_labeller(nlabs))) +
    ggtitle("Power vs. CV, ICC = 0.01") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x  = element_text(angle = 90)) +
    guides(color = FALSE))

(dfplot2 <- ggplot(df2, aes(x = cv, y = power)) + ylim(0,1) +
    geom_abline(aes(slope = 0, intercept = 0.8, color = "red")) +
    geom_linerange(aes(ymin = low, ymax = high), size = 3) +
    facet_grid(n ~ m, labeller = labeller(m = as_labeller(mlabs),
                                          n = as_labeller(nlabs))) +
    ggtitle("Power vs. CV, ICC = 0.05") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x  = element_text(angle = 90)) +
    guides(color = FALSE))
