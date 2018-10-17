dev.off()
rm(list = ls())
cat("\014")
set.seed(56)

# ============================================= #
# ===== 4 ways of simulating GBM ============== #
# ============================================= #

# ===== Global parameters ===================== #
spot       <- 100
rate       <- 0.05
sig        <- 0.2
mat        <- 10
time_steps <- 52
sims       <- 100000
dt         <- 1 / time_steps
periods    <- mat * time_steps

# Split sims in two half for antithetic variates
half_sim   <- sims/2
rest_sim   <- half_sim + 1
# ============================================= #

# ============================================= #
# ===== 1. Loop with antithetic variates ====== #
# ============================================= #

ptm <- proc.time()
z <- rnorm(sims * periods / 2, 0, 1)
z <- matrix(z, sims / 2, periods)

drift <- (rate - 0.5 * sig ** 2) * dt
vol   <- sig * sqrt(dt)

stock1 <- matrix(spot, sims, periods)
for (t in 2:(periods)){
  stock1[1:half_sim,t] <- 
    stock1[1:half_sim,t-1] * exp(drift + vol * z[,t])
  stock1[rest_sim:sims,t] <- 
    stock1[rest_sim:sims,t-1] * exp(drift - vol * z[,t])
}
time1 <- proc.time() - ptm
# ============================================= #

# ============================================= #
# ===== 2. Loop with predefined exponentials == #
# ============================================= #
ptm   <- proc.time()
z <- rnorm(sims * periods / 2, 0, 1)
z <- matrix(z, sims / 2, periods)

drift <- (rate - 0.5 * sig ** 2) * dt
vol   <- sig * sqrt(dt) * z

exp1  <- exp(drift + vol)
exp2  <- exp(drift - vol)

stock2 <- matrix(spot, sims, periods)
for (t in 2:(periods)){
  stock2[1:half_sim,t]    <- stock2[1:half_sim,t-1] * exp1[,t]
  stock2[rest_sim:sims,t] <- stock2[rest_sim:sims,t-1] * exp2[,t]
}
time2 <- proc.time() - ptm
# ============================================= #

# ============================================= #
# ===== 3. Remove for loop using apply fct ==== #
# ============================================= #
ptm   <- proc.time()
z <- rnorm(sims * periods / 2, 0, 1)
z <- matrix(z, sims / 2, periods)

drift <- (rate - 0.5 * sig ** 2) * dt
vol   <- sig * sqrt(dt)

stock3 <- matrix(log(spot), sims, periods)
stock3[1:half_sim, 2:periods] <- 
  log(spot) + t(apply(drift + vol * z[,2:periods], 1, cumsum))
stock3[rest_sim:sims, 2:periods] <- 
  log(spot) + t(apply(drift - vol * z[,2:periods], 1, cumsum))
stock3 <- exp(stock3)
time3 <- proc.time() - ptm
# ============================================= #

# ============================================= #
# ===== 4. For loop, no antithetic variates === #
# ============================================= #
ptm <- proc.time()
z <- rnorm(sims * periods, 0, 1)
z <- matrix(z, sims, periods)

drift <- (rate - 0.5 * sig ** 2) * dt
vol   <- sig * sqrt(dt)

stock4 <- matrix(spot, sims, periods)
for (t in 2:(periods)){
  stock4[,t] <- stock4[,t-1] * exp(drift + vol * z[,t])
}
time4 <- proc.time() - ptm
# ============================================= #

# ============================================= #
# ===== No number generation comparisons ====== #
# ============================================= #

z <- rnorm(sims * periods / 2, 0, 1)
z <- matrix(z, sims / 2, periods)

# ============================================= #
# ===== 5. Loop with antithetic variates ====== #
# ============================================= #

ptm <- proc.time()

drift <- (rate - 0.5 * sig ** 2) * dt
vol   <- sig * sqrt(dt)

stock1 <- matrix(spot, sims, periods)
for (t in 2:(periods)){
  stock1[1:half_sim,t] <- 
    stock1[1:half_sim,t-1] * exp(drift + vol * z[,t])
  stock1[rest_sim:sims,t] <- 
    stock1[rest_sim:sims,t-1] * exp(drift - vol * z[,t])
}
time5 <- proc.time() - ptm
# ============================================= #

# ============================================= #
# ===== 6. Loop with predefined exponentials == #
# ============================================= #
ptm   <- proc.time()

drift <- (rate - 0.5 * sig ** 2) * dt
vol   <- sig * sqrt(dt) * z

exp1  <- exp(drift + vol)
exp2  <- exp(drift - vol)

stock2 <- matrix(spot, sims, periods)
for (t in 2:(periods)){
  stock2[1:half_sim,t]    <- stock2[1:half_sim,t-1] * exp1[,t]
  stock2[rest_sim:sims,t] <- stock2[rest_sim:sims,t-1] * exp2[,t]
}
time6 <- proc.time() - ptm
# ============================================= #

# ============================================= #
# ===== 7. Remove for loop using apply fct ==== #
# ============================================= #
ptm   <- proc.time()

drift <- (rate - 0.5 * sig ** 2) * dt
vol   <- sig * sqrt(dt)

stock3 <- matrix(log(spot), sims, periods)
stock3[1:half_sim, 2:periods] <- 
  log(spot) + t(apply(drift + vol * z[,2:periods], 1, cumsum))
stock3[rest_sim:sims, 2:periods] <- 
  log(spot) + t(apply(drift - vol * z[,2:periods], 1, cumsum))
stock3 <- exp(stock3)
time7 <- proc.time() - ptm
# ============================================= #

# ============================================= #
# ===== Conclusions =========================== #
# ============================================= #

# Considering random number generation
time1
time2
time3
time4

# Without considering random number generation
time5
time6
time7

# ============================================= #
#' When considering random number generation
#' the first for loop is the fastest.
#' When not considering random number generation
#' the predefining method becomes faster.
#' This is a bit strange. Since random numbers
#' do need generating, I will rely on the first
#' set of results.
# ============================================= #