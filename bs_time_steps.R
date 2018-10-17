dev.off()
rm(list = ls())
cat("\014")
options(scipen = 999)
library(pbivnorm)
source("C:/Users/wigr11ab/Dropbox/KU/K3/MCM/bs_function.R")

# ============================================= #
# ===== Global parameters ===================== #
# ============================================= #

sims    <- 100000
mat     <- 10
steps   <- 252
periods <- mat * steps
dt      <- 1 / steps

spot   <- 100
r      <- 0.05
sig    <- 0.2
strike <- 100

# ============================================= #
# ===== Building the wiener process frame ===== #
# ============================================= #

av_mc <- function(spot = 100, r = 0.05, sig = 0.2, mat = 10, steps = 1, strike = 100, sims = 50000){
  periods <- steps * mat
  dt      <- 1 / steps
  
  z <- rnorm(sims * periods)
  z <- matrix(z, sims, periods)
  
  w  <- rep(spot, sims * periods)
  w1 <- matrix(w, sims, periods)
  w2 <- matrix(w, sims, periods)
  
  drift <- (r - 0.5 * sig ** 2) * dt
  vol   <- sig * sqrt(dt)
  
  for(t in 2:(periods)){
    w1[, t] <- w1[, t-1] * exp(drift + vol * z[,t])
    w2[, t] <- w2[, t-1] * exp(drift - vol * z[,t])
  }
  
  dc     <- exp( -r * mat)
  p_off1 <- dc * pmax(w1[,periods] - strike, 0)
  p_off2 <- dc * pmax(w2[,periods] - strike, 0)
  p_off  <- 0.5 * (p_off1 + p_off2)
  
  price <- sum(p_off) / sims
  sd    <- sqrt(sum((p_off - price) ** 2) / (sims - 1))
  se    <- 1.96 * sd / sqrt(sims)
  lower <- price - se
  upper <- price + se
  
  return(c(price = price, se = se, lower = lower, upper = upper))
}

av_mc2 <- function(spot = 100, r = 0.05, sig = 0.2, mat = 10, strike = 100, sims = 50000){
  z <- rnorm(sims)

  drift <- (r - 0.5 * sig ** 2) * mat
  vol   <- sig * sqrt(mat)
  
  for(t in 2:(periods)){
    w1 <- spot * exp(drift + vol * z)
    w2 <- spot * exp(drift - vol * z)
  }
  
  dc    <- exp( -r * mat)
  
  p_off1 <- dc * pmax(w1 - strike, 0)
  p_off2 <- dc * pmax(w2 - strike, 0)
  p_off  <- 0.5 * (p_off1 + p_off2)

  price <- sum(p_off) / sims
  sd    <- sqrt(sum((p_off - price) ** 2) / (sims - 1))
  se    <- 1.96 * sd / sqrt(sims)
  lower <- price - se
  upper <- price + se
  
  return(c(price = price, se = se, lower = lower, upper = upper))
}

av_mc(steps = steps)
av_mc2()
bs(price = spot, rate = r, vol = sig, mat = 10, strike = 100)

w <- rep(spot, sims * periods)
w <- matrix(w, sims, periods)

# ===== Building the standard normal randoms == #
z <- rnorm(sims / 2 * periods, 0, 1)
z <- matrix(z, sims / 2, periods)
# z <- rnorm(sims * periods * steps, 0, 1)
# z <- matrix(z, sims, periods * steps)


# ===== Building the simulated wiener paths === #
drift <- (r - 0.5 * v ** 2) * dt
vol   <- v * sqrt(dt)

# Implementing antithetic variates
plus    <- c(1:(sims / 2))
anti    <- c((1 + sims / 2):sims)

test1   <- spot * exp(drift + vol * z[,periods])
test2   <- spot * exp(drift - vol * z[,periods])
test    <- 0.5 * (test1 + test2)

payoff <- dc * pmax(test - strike, 0)
mean_p <- mean(payoff)
sd_p   <- sd(payoff)
se_p   <- 1.96 * sd / sqrt(sims / 2)
lower  <- mean_p - se_p
upper  <- mean_p + se_p

c(mean_p, lower, upper, se_p)

for(t in 2:(periods)){
  w[plus, t] <- w[plus, t-1] * exp(drift + vol * z[,t])
  w[anti, t] <- w[anti, t-1] * exp(drift - vol * z[,t])
  # w[, t] <- w[, t-1] * exp(drift + vol * z[,t])
}
w <- 0.5 * (w[plus,] + w[anti,])

# ===== Finding the BS-price ================== #
s_T <- w[, periods]

c_T <- pmax(s_T - strike, 0)
c_h <- sum(c_T) / (sims / 2)

dc   <- exp(-r * mat)
c_t  <- dc * c_T
c_ht <- dc * c_h

# ===== Deriving the interval ================= #
sd    <- sqrt(1 / (sims/2 - 1) * sum((c_t - c_ht) ** 2))
se    <- 1.96 * sd / sqrt(sims / 2)
lower <- c_ht - se
upper <- c_ht + se

interval_cmc <- c(lower = lower, 
                  value = c_ht,
                  upper = upper,
                  se    = se)

# ===== Control variate ================================= #
beta  <- sum((c_T - c_h) * (s_T - mean(s_T))) / sum((s_T - mean(s_T)) ** 2)
c_0   <- dc * sum(c_T - beta * (s_T - spot * exp(r * mat))) / (sims / 2)
sb    <- sqrt(var(c_T) / (sims / 2 - 1)) # only works asymptotically - is it right?
lower <- c_0 - 1.96 * sb / sqrt(sims / 2)
upper <- c_0 + 1.96 * sb / sqrt(sims / 2)

interval_cv <- c(lower = lower,
                 value = c_0,
                 upper = upper,
                 se    = sb / sqrt(sims))

# Compare crude monte carlo (cmc) with control variate (cv)
interval_cmc
interval_cv

# ===== The analytical option price =========== #
bs(price = spot, vol = v, mat = mat, rate = r, strike = strike)

#' Conclusion:
#' It seems my control variate (cv) estimate misses the target, even
#' if it is much closer than the crude monte carlo (cmc) method. The
#' interval does not contain the true price.
#' 
#' I am also not seeing any reduction in variance from using
#' antithetic variates.