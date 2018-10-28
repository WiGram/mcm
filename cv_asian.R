dev.off()             
rm(list = ls())       
cat("\014")           
options(scipen = 999) 
library(stats4)       
library(ggplot2)      
set.seed(12345)

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

n <- 100000
s <- 100
strike <- seq(90, 110, 10)
# rate   <- seq(0.03,0.07,0.02)
# strike <- 55
rate   <- 0.05
# sigma  <- seq(0.2,0.40,0.1)
sigma <- 0.2

mat <- 1        # e.g. four months
ts  <- 500       # time steps, e.g. 22 per month
dt  <- mat / ts

z <- rnorm(ts * n, 0, 1)
z <- matrix(z, nrow = n, ncol = ts)
S <- matrix(s, nrow = n, ncol = ts)

prices_cmc <- sds_cmc <- list()
prices_cv  <- sds_cv   <- list()

price_cmc <- sd_cmc <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_cv  <- sd_cv  <- matrix(0, nrow = length(strike), ncol = length(sigma))

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
      c_g_bar <- c_g_bar_fct(s, k, r, v, mat)
      
      beta <- beta_fct(c_a, c_a_bar, c_g, c_g_bar)
      
      price_cmc[h, i] <- c_a_bar
      sd_cmc[h, i] <- sd(c_a) / sqrt(n)
      
      price_cv[h, i] <- mean(c_a - beta * (c_g - c_g_bar))
      sd_cv[h, i]    <- sd(c_a - beta * (c_g - c_g_bar)) / sqrt(n)
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
  
  prices_cv[[a]] <- price_cv
  rownames(prices_cv[[a]]) <- paste0('k: ', strike)
  colnames(prices_cv[[a]]) <- paste0('v: ', sigma)
  
  sds_cv[[a]] <- sd_cv
  rownames(sds_cv[[a]]) <- paste0('k: ', strike)
  colnames(sds_cv[[a]]) <- paste0('v: ', sigma)
}

prices_cmc
prices_cv