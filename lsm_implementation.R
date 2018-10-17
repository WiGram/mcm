source('C:/Users/wigr11ab/Dropbox/KU/K3/MCM/lsm_function.R')
# lsm <- function(sim_matrix, strike, spot, r)

source('C:/Users/wigr11ab/Dropbox/KU/K3/MCM/gbm_function.R')
# gbm <- function(spot, interest, volatility, paths, periods, anti = TRUE, seed = 0)

# EXAMPLE PARAMETERS
strike  <- 1.1
spot    <- 1.0
r       <- 0.06

EXAMPLE = matrix( c(1, 1.09, 1.08, 1.34,
                    1, 1.16, 1.26, 1.54,
                    1, 1.22, 1.07, 1.03,
                    1, 0.93, 0.97, 0.92,
                    1, 1.11, 1.56, 1.52,
                    1, 0.76, 0.77, 0.90,
                    1, 0.92, 0.84, 1.01,
                    1, 0.88, 1.22, 1.34), 
                  nrow=8, ncol=4, byrow=TRUE)

# Simulation parameters
vol     <- 0.20
paths   <- 8
periods <- 4

# The LSM function applies the LSM algorithm and returns the price.
lsm(EXAMPLE, strike, spot, r)

# The same as before, but with simulated matrix
sim_anti_T <- gbm(spot, r, vol, 100, 10, anti = TRUE, seed = 12345)
sim_anti_F <- gbm(spot, r, vol, 100, 10, anti = FALSE, seed = 12345)

lsm(sim_anti_T, strike, spot, r)
lsm(sim_anti_F, strike, spot, r)


# Finally repeating the simulated experiment 100 times
reps <- 100

lsm_rep_aT <- list()
for (i in 1:reps){
  sim <- gbm(spot, r, vol, 100, 10, anti = TRUE)
  lsm_rep_aT[[i]] <- lsm(sim, strike, spot, r)
}
lsm_rep_aT <- do.call(rbind, lsm_rep_aT)
price_avg_aT <- mean(lsm_rep_aT)
price_std_aT <- sd(lsm_rep_aT)

bounds_T <- c(price_avg_aT - 1.96 * price_std_aT / sqrt(reps), price_avg_aT, price_avg_aT + 1.96 * price_std_aT / sqrt(reps))

lsm_rep_aF <- list()
for (i in 1:reps){
  sim <- gbm(spot, r, vol, 100, 10, anti = FALSE)
  lsm_rep_aF[[i]] <- lsm(sim, strike, spot, r)
}
lsm_rep_aF <- do.call(rbind, lsm_rep_aF)
price_avg_aF <- mean(lsm_rep_aF)
price_std_aF <- sd(lsm_rep_aF)

bounds_F <- c(price_avg_aF - 1.96 * price_std_aF / sqrt(reps), price_avg_aF, price_avg_aF + 1.96 * price_std_aF / sqrt(reps))
