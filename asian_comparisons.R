dev.off()             
rm(list = ls())       
cat("\014")           
options(scipen = 999) 
library(stats4)       
library(ggplot2)      

# ============================================= #
# ===== Asian options ========================= #
# ============================================= #

n <- 10000
s <- 100
strike <- seq(90, 110, 10)
rate   <- seq(0.03,0.07,0.02)
sigma  <- seq(0.2,0.40,0.1)

mat <- 1        # e.g. four months
ts  <- 250       # time steps, e.g. 22 per month
dt  <- mat / ts

z <- rnorm(ts * n, 0, 1)
z <- matrix(z, nrow = n, ncol = ts)
S <- matrix(s, nrow = n, ncol = ts)

prices_a <- sds_a <- list()
prices_g <- sds_g <- list()

price_a <- sd_a <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_g <- sd_g <- matrix(0, nrow = length(strike), ncol = length(sigma))

h <- 1
i <- 1

for (r in rate){
  for (v in sigma){
    for (k in strike){
      
      drift <- (r - 0.5 * v ** 2) * dt
      vol   <- v * sqrt(dt)
      
      for (t in 2:ts){
        S[,t] <- S[, t-1] * exp(drift + vol * z[,t])
      }
      
      arith_mean <- rowMeans(S)
      geom_mean  <- exp(rowMeans(log(S)))
      
      c_a <- exp(-r * mat) * pmax(arith_mean - k, 0)
      c_g <- exp(-r * mat) * pmax(geom_mean  - k, 0)
      
      c_a_bar <- mean(c_a)
      c_g_bar <- mean(c_g)
      
      price_a[h, i] <- c_a_bar
      sd_a[h, i] <- sd(c_a) / sqrt(n)
      
      price_g[h, i] <- c_g_bar
      sd_g[h, i] <- sd(c_g) / sqrt(n)
      
      h <- h + 1
    }
    h <- 1
    i <- i + 1
  }
  h <- 1
  i <- 1
  a <- paste0('r: ', r)
  prices_a[[a]] <- price_a
  rownames(prices_a[[a]]) <- paste0('k: ', strike)
  colnames(prices_a[[a]]) <- paste0('v: ', sigma)
  
  prices_g[[a]] <- price_g
  rownames(prices_g[[a]]) <- paste0('k: ', strike)
  colnames(prices_g[[a]]) <- paste0('v: ', sigma)
  
  sds_a[[a]] <- sd_a
  rownames(sds_a[[a]]) <- paste0('k: ', strike)
  colnames(sds_a[[a]]) <- paste0('v: ', sigma)

  sds_g[[a]] <- sd_g
  rownames(sds_g[[a]]) <- paste0('k: ', strike)
  colnames(sds_g[[a]]) <- paste0('v: ', sigma)
}