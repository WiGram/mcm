dev.off()             
rm(list = ls())       
cat("\014")           
options(scipen = 999) 
library(stats4)       
library(ggplot2)      

# ============================================= #
# ===== Asian options ========================= #
# ============================================= #

# ----- Functions ----------------------------- #
# Kemma and Vorst
c_g_bar_fct <- function(s, k, r, v, mat){
  d_star <- 0.5 * (r - v ** 2 / 6) * mat
  d1     <- (log(s / k) + 0.5 * (r + v ** 2 / 6) * mat) / (v * sqrt(mat / 3))
  d2     <- d1 - v * sqrt(mat / 3)
  exp(-r * mat) * (exp(d_star) * s * pnorm(d1) - k * pnorm(d2))
}

# CV beta
beta_fct <- function(c_1, c_1_bar, c_2, c_2_bar){
  sum((c_1 - c_1_bar) * (c_2 - c_2_bar)) /
    sum((c_2 - c_2_bar) ** 2)
}
# --------------------------------------------- #

set.seed(12345)

n <- 10000
s <- 100
strike <- seq(90, 110, 10)
rate   <- seq(0.03,0.07,0.02)
# strike <- 55
# rate   <- 0.05
sigma  <- seq(0.2,0.40,0.1)
# sigma <- 0.2

mat <- 1        # e.g. four months
ts  <- 250       # time steps, e.g. 22 per month
dt  <- mat / ts

z <- rnorm(ts * n, 0, 1)
z_av <- z[1:(length(z) / 2)]

z    <- matrix(z, nrow = n, ncol = ts)
z_av <- matrix(z, nrow = n / 2, ncol = ts)

S  <- matrix(s, nrow = n, ncol = ts)
S1 <- matrix(s, nrow = n / 2, ncol = ts)
S2 <- matrix(s, nrow = n / 2, ncol = ts)

prices_cmc  <- sds_cmc   <- list()
prices_av   <- sds_av    <- list()
prices_cv   <- sds_cv    <- list()
prices_cav  <- sds_cav   <- list()

price_cmc  <- sd_cmc  <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_av   <- sd_av   <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_cv   <- sd_cv   <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_cav  <- sd_cav  <- matrix(0, nrow = length(strike), ncol = length(sigma))


h <- 1
i <- 1

for (r in rate){
  for (v in sigma){
    for (k in strike){
      
      drift <- (r - 0.5 * v ** 2) * dt
      vol   <- v * sqrt(dt)
      
      for (t in 2:ts){
        S[, t]  <- S[, t-1] * exp(drift + vol * z[, t])
        S1[,t] <- S1[, t-1] * exp(drift + vol * z_av[,t])
        S2[,t] <- S2[, t-1] * exp(drift - vol * z_av[,t])
      }
      
      # CMC
      mean_a   <- rowMeans(S)
      
      c_a      <- exp(-r * mat) * pmax(mean_a - k, 0)
      c_a_bar  <- mean(c_a)

      price_cmc[h, i] <- c_a_bar
      sd_cmc[h, i]    <- sd(c_a) / sqrt(n)
            
      #AV
      mean_a_1 <- rowMeans(S1)
      mean_a_2 <- rowMeans(S2)
      
      c_a_1    <- exp(-r * mat) * pmax(mean_a_1 - k, 0)
      c_a_2    <- exp(-r * mat) * pmax(mean_a_2 - k, 0)
      c_av     <- 0.5 * (c_a_1 + c_a_2)
      c_av_bar <- mean(c_av)
      
      price_av[h, i] <- c_av_bar
      sd_av[h, i]    <- sd(c_av) / sqrt(n)
      
      # CV
      mean_g  <- exp(rowMeans(log(S)))
      
      c_g     <- exp(-r * mat) * pmax(mean_g - k, 0)
      c_g_bar <- c_g_bar_fct(s, k, r, v, mat)
      beta_cv <- beta_fct(c_a, c_a_bar, c_g,  c_g_bar)
      
      price_cv[h, i]  <- mean(c_a - beta_cv * (c_g - c_g_bar))
      sd_cv[h, i]     <- sd(c_a - beta_cv * (c_g - c_g_bar)) / sqrt(n)
      
      # CV + AV (_ca)
      mean_g_1  <- exp(rowMeans(log(S1)))
      mean_g_2  <- exp(rowMeans(log(S2)))
      
      c_g_1   <- exp(-r * mat) * pmax(mean_g_1  - k, 0)
      c_g_2   <- exp(-r * mat) * pmax(mean_g_2  - k, 0)
      c_ga    <- 0.5 * (c_g_1 + c_g_2)
      c_g_bar <- c_g_bar
      beta_ca <- beta_fct(c_av, c_av_bar, c_ga, c_g_bar)
    
      price_cav[h, i] <- mean(c_av - beta_ca * (c_ga - c_g_bar))
      sd_cav[h, i]    <- sd(c_av - beta_ca * (c_ga - c_g_bar)) / sqrt(n)
      h <- h + 1
    }
    h <- 1
    i <- i + 1
  }
  h <- 1
  i <- 1
  a <- paste0('r: ', r)
  prices_cmc[[a]] <- price_cmc
  rownames(prices_cmc[[a]]) <- paste0('k: ', strike)
  colnames(prices_cmc[[a]]) <- paste0('v: ', sigma)
  
  sds_cmc[[a]] <- sd_cmc
  rownames(sds_cmc[[a]]) <- paste0('k: ', strike)
  colnames(sds_cmc[[a]]) <- paste0('v: ', sigma)
  
  prices_av[[a]] <- price_av
  rownames(prices_av[[a]]) <- paste0('k: ', strike)
  colnames(prices_av[[a]]) <- paste0('v: ', sigma)
  
  sds_av[[a]] <- sd_av
  rownames(sds_av[[a]]) <- paste0('k: ', strike)
  colnames(sds_av[[a]]) <- paste0('v: ', sigma)
  
  prices_cv[[a]] <- price_cv
  rownames(prices_cv[[a]]) <- paste0('k: ', strike)
  colnames(prices_cv[[a]]) <- paste0('v: ', sigma)
  
  sds_cv[[a]] <- sd_cv
  rownames(sds_cv[[a]]) <- paste0('k: ', strike)
  colnames(sds_cv[[a]]) <- paste0('v: ', sigma)

  prices_cav[[a]] <- price_cav
  rownames(prices_cav[[a]]) <- paste0('k: ', strike)
  colnames(prices_cav[[a]]) <- paste0('v: ', sigma)
  
  sds_cav[[a]] <- sd_cav
  rownames(sds_cav[[a]]) <- paste0('k: ', strike)
  colnames(sds_cav[[a]]) <- paste0('v: ', sigma)
}