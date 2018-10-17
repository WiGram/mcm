dev.off()
rm(list = ls())
cat("\014")
library(pbivnorm)

# ================================================= #
# =========== Pricing a compound option =========== #
# ================================================= #

source("C:/Users/wigr11ab/Dropbox/KU/K3/MCM/bs_function.R")

# Simulate n paths for T1
# Simulate k paths for T2 for each
# stock path simulated at T1
# We need only the final period stock price
s    <- 100
n    <- 100
m    <- 1000000
r    <- 0.05
v    <- 0.2
k_1  <- 10
k_2  <- 100
mat1 <- 1
mat2 <- 3

# solve_problem, optim minimises, unles control = list(fnscale = -1)


compound_bs <- function(spot, vol, rate, mat1, mat2, strike1, strike2){
  
  # First we need to find the tipping price, bar(s)
  optim_s_bar <- function(theta, vol, rate, mat, strike1, strike2){
    s_bar <- theta
    
    k2 <- log(s_bar / strike1) + (rate - 0.5 * vol ** 2) * mat
    k2 <- k2 / (vol * sqrt(mat))
    k1 <- k2 + vol * sqrt(mat)
    
    dc    <- exp(-rate * mat)
    n1    <- pnorm(k1)
    n2    <- pnorm(k2)
    
    min   <- abs(s_bar * n1 - strike1 * dc * n2 - strike2)
    
    return(min)
  }
  
  # By setting maximum = FALSE the solver minimises
  s_bar <- optimise(f = optim_s_bar, 
                    interval = c(0, 200),
                    maximum = FALSE,
                    vol = vol, 
                    rate = rate, 
                    mat = mat1, 
                    strike1 = strike1, 
                    strike2 = strike2)
  
  s_bar <- s_bar$minimum
  
  h2 <- log(spot / s_bar) + (rate - 0.5 * vol ** 2) * mat1
  h2 <- h2 / (vol * sqrt(mat1))
  h1 <- h2 + vol * sqrt(mat1)
  
  k2 <- log(spot / strike2) + (rate - 0.5 * vol ** 2) * mat2
  k2 <- k2 / (vol * sqrt(mat2))
  k1 <- k2 + vol * sqrt(mat2)
  
  corr <- sqrt(mat1 / mat2)
  
  n3 <- pbivnorm(h1, k1, corr)
  n2 <- pbivnorm(h2, k2, corr)
  n1 <- pnorm(h2)
  
  dc2 <- exp(-rate * mat2)
  dc1 <- exp(-rate * mat1)
  
  c <- spot * n3 - strike2 * dc2 * n2 - strike1 * dc1 * n1
  
  return(c)
}

compound_bs(spot = s, vol = v, rate = r, 
            mat1 = mat1, mat2 = mat2,
            strike1 = k_1, strike2 = k_2)

compound_mc <- function(sim1, sim2, s, r, v, k_1, k_2, mat1, mat2){
  # sim1 = n
  # sim2 = m
  
  # Discount factors
  dc_1  <- exp(-r * mat1)
  dc_2  <- exp(-r *  (mat2 - mat1))
  
  # Price matrices
  s_1 <- rep(s, sim1)
  s_2 <- matrix(0,   nrow = sim1, ncol = sim2)
  
  # Random numbers
  z   <- rnorm(n = sim1, 0, 1)
  z_2 <- rnorm(n = sim1 * sim2, 0, 1)
  z_2 <- matrix(z_2, nrow = sim1, ncol = sim2)
  
  # First, compute T1 stock prices for n stocks
  drift <- (r - 0.5 * v ** 2) * mat1
  vol   <- v * sqrt(mat1)
  s_1   <- s_1 * exp(drift + vol * z)
  
  # For each n stock price at T1, m stock prices are computed for T2
  drift <- (r - 0.5 * v ** 2) * (mat2 - mat1)
  vol   <- v * sqrt(mat2 - mat1)
  for (i in 1:sim1){
    s_2[i,] <- s_1[i] * exp(drift + vol * z_2[i,])
  }
  
  # The second option price is the T1-PV of the average of each of m paths
  c_2 <- 1 / sim2 * dc_2 * rowSums(pmax(s_2 - k_2, 0))
  
  # The final option price is the T0-PV of the average of each of n paths
  c_1 <- dc_1 * pmax(c_2 - k_1, 0)
  C   <- sum(c_1) / sim1
  
  se  <- sqrt(1 / (sim1 - 1) * sum((c_1 - C) ** 2))
  se  <- se / sqrt(sim1)
  
  lower = round(C - 1.96 * se, 4)
  upper = round(C + 1.96 * se, 4)
  
  int <- c(lower = lower, 
           estimate = round(C,4), 
           upper = upper,
           se = round(se, 4))
  
  # ================================================= #
  # Applying BS-pricing formula to compare with bias.
  # ================================================= #
  C_1 <- compound_bs(spot = s, vol = v, rate = r, 
                     mat1 = mat1, mat2 = mat2,
                     strike1 = k_1, strike2 = k_2)
  
  abs_bias <- round(C - C_1, 4)
  rel_bias <- round((C - C_1) / C_1 * 100, 4)
  
  return(list(MC_price = int,
              BS_price = round(C_1,4),
              abs_bias = abs_bias,
              rel_bias = rel_bias))
}

prices <- compound_mc(n, m, s, r, v, k_1, k_2, mat1, mat2)

#' Conclusion: We have now priced a compound option.
#' The compound option is biased high in this example.
#' Using the BS-formula to derive analytical solutions, 
#' we confirm the MC method is bias high.

paste0("MC Standard error is: ", prices$MC_price[4], ".", sep = "")

paste0("Black Scholes price is: ", prices$BS_price, " USD.", sep = "")

paste0("MC price is: ", prices$MC_price[2], " USD.", sep = "")

paste0("Surcharge of: ", prices$abs_bias, " USD on actual price.", sep = "")

# Or in percent
paste0("Surcharge of: ", prices$rel_bias, " percent of actual price.", sep = "")
# ================================================= #
