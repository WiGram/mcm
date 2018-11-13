# Pricing European options using Black Scholes and MC method

# The following source contains the black scholes formula
# (Not needed - I wrote the function into this script directly)
# source("C:/Users/wigr11ab/Dropbox/KU/K3/MCM/bs_function.R")

library(ggplot2)

# ================================== #
# ===== Black Scholes method ======= #
# ================================== #

# Black Scholes parameters
price  <- 100
vol    <- 0.2
mat    <- 1
rate   <- 0.05
strike <- 110

# Black Scholes option pricing formula
bs <- function(price, vol, mat, rate, strike){
  d1 <- (log(price / strike) + (rate + 0.5 * vol ** 2) * mat)/(vol * sqrt(mat))
  d2 <- d1 - vol * sqrt(mat)
  
  # mean and sd do not need to be specified: default to N(0,1)
  n1 <- pnorm(d1, mean = 0, sd = 1)
  n2 <- pnorm(d2, mean = 0, sd = 1)
  dc <- exp(-rate * mat)
  
  return(price * n1 - dc * strike * n2)
}

# ================================== #
# ===== Monte Carlo method ========= #
# ================================== #

# Random number parameters
sims   <- 1000
z_mean <- 0
z_sd   <- 1

# Monte carlo function
bs_mc <- function(sims, 
                  z_mean, z_sd,
                  price, vol, mat, rate, strike)
  {
  z  <- rnorm(n = sims, mean = z_mean, sd = z_sd)
  
  s0 <- price
  s  <- s0 * exp((rate - 0.5 * vol ** 2) * mat + vol * sqrt(mat) * z)
  
  dc  <- exp(-rate * mat)
  c   <- dc * pmax( s - strike, 0)
  c_T <- sum(c) / sims
  
  se  <- sqrt(1 / (sims - 1) * sum((c - c_T) ** 2))
  
  int <- c(lower = c_T - 1.96 * se / sqrt(sims),
           estimate = c_T,
           upper = c_T + 1.96 * se / sqrt(sims))
  
  return(list(stock = s,
              option = c,
              estimates = int, 
              se = se))
}

# Monte carlo parameter estimation and interval
mc_sim <- bs_mc(sims, 
                z_mean, z_sd, 
                price, vol, mat, rate, strike)
mc_sim[3:4]

ggplot(NULL) +
  geom_point(aes(x = mc_sim[[1]],
                y = mc_sim[[2]]),
             shape = 1, alpha = 0.2,
             col = 'blue', size = 2) +
  labs(x = 'Stock price', y = 'Option price')

# Black Scholes price
bs <- bs(price, vol, mat, rate, strike)
bs



# ============================== #
# ===== Path dependence (pd) === #
# ============================== #

# Final parameters
periods    <- 10
time_steps <- 10

#' Path dependent MC simulation e.g. Asian options
#' where Asian option payoffs depend on average
#' level of the underlying asset over some period
pd_mc <- function(sims, periods, time_steps,
                  z_mean, z_sd,
                  price, vol, mat, rate, strike)
  {
  s0 <- price
  dt <- 1 / time_steps
  
  z  <- rnorm(n = sims * periods * time_steps, mean = 0, sd = 1)
  z  <- matrix(z, nrow = sims, ncol = periods * time_steps)
  
  s  <- matrix(s0, nrow = sims, ncol = periods * time_steps)
  
  for (t in 2:(periods * time_steps)){
    s[,t] <- s[,t-1] * exp((rate - 0.5 * vol ** 2) * dt + 
                             vol * sqrt(dt) * z[,t])
  }
  
  s_bar <- rowMeans(s)
  c_bar <- exp(-rate * mat) * pmax(s_bar - strike, 0)
  
  c_hat <- mean(c_bar)
  c_hat
  
  se  <- sqrt(1 / (sims - 1) * sum((c_bar - c_hat) ** 2))
  se  <- se / sqrt(sims)
  
  int <- c(lower = c_hat - 1.96 * se,
           estimate = c_hat,
           upper = c_hat + 1.96 * se)
  
  return(list(estimates = int,
              se = se))
}

pd_mc(sims, periods, time_steps,
      z_mean, z_sd,
      price, vol, mat, rate, strike)