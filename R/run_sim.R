# Run a simulation

# This file contains the code needed to generate the data for the paper,
#  as well as generating the power versus coefficient of variation plots.

# To run this simulation insure that the packages below are installed, and
#  insure that your working directory contains this file, run_sim.R, as well
#  as the files sourced below.

library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(latex2exp)

setwd("~/crtpower/R")

source("misc_functions.R") # miscellaneous functions
source("crt_power_sim.R") # function to calculate power by simulation
source("power_crt_test.R") # function for power analysis calculations
source("get_df.R") # function to create data frame to iterate over
source("sim_power.R") # function for running the simulation

# initialize vectors to pass to get_df() function
ICC0 <- c(0.01, 0.05)
m0 <- c(5, 25, 50)
n0 <- c(10, 30, 100)
d0 <- 0
cv0 <- seq(0,0.55, by = 0.05)

# create data frame to iterate over
df <- get_df(ICC0, m0, n0, cv0, d0)

nsim <- 20000

# string used when warn = TRUE in sim_power() below
warning_string <- "bobyqa"

# counter for keeping track of bobyqa warnings
ww <- numeric(nrow(df))

set.seed(2001)

# simulate power
power_list <- sim_power(nsim,
                        m = df$m, n = df$n, a = df$a,
                        d = df$d, ICC = df$ICC,
                        warn = TRUE)

# store the output of the simulation in df
df$low <- power_list$low
df$power <- power_list$power
df$high <- power_list$high

# labels for m and n
mlabs <- c(`5` = "M = 10", `25` = "M = 50", `50` = "M = 100")
nlabs <- c(`10` = TeX("$\\bar{n} = 10$"),
           `30` = TeX("$\\bar{n} = 30$"),
           `100` = TeX("$\\bar{n} = 100$"))

df1 <- read.csv("crtsim_20000.csv")


# subset by ICC
df1 <- filter(df, ICC == 0.01)
df2 <- filter(df, ICC == 0.05)

# generate plots
(dfplot1 <- ggplot(df1, aes(x = cv, y = power)) + ylim(0,1) +
    geom_abline(aes(slope = 0, intercept = 0.8, color = "red")) +
    geom_linerange(aes(ymin = low, ymax = high), size = 3) +
    facet_grid(n ~ m,
               labeller = labeller(m = as_labeller(mlabs),
                                   n = as_labeller(nlabs, default = label_parsed))) +
    ggtitle("Power vs. CV, ICC = 0.01") +
    xlab("Coefficient of Variation") + ylab("Power") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x  = element_text(angle = 90)) +
    guides(color = FALSE))

ggsave("crtsim_20000_1.pdf", plot = dfplot1, height = 10, width = 7)

(dfplot2 <- ggplot(df2, aes(x = cv, y = power)) + ylim(0,1) +
    geom_abline(aes(slope = 0, intercept = 0.8, color = "red")) +
    geom_linerange(aes(ymin = low, ymax = high), size = 3) +
    facet_grid(n ~ m, labeller = labeller(m = as_labeller(mlabs),
                                          n = as_labeller(nlabs,
                                                          default = label_parsed)
    )) +
    ggtitle("Power vs. CV, ICC = 0.05") +
    xlab("Coefficient of Variation") + ylab("Power") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x  = element_text(angle = 90)) +
    guides(color = FALSE))

ggsave("crtsim_20000_5.pdf", plot = dfplot2, height = 10, width = 7)