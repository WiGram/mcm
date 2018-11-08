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

# Brownian bridge
bridge <- function(data, stop, n){
  for (t in 2:(stop-1)){
    a <- data[,(t-1)]
    b <- data[,stop]
    data[,t] <- ((stop - t) * a + b) / (stop - (t-1)) +
      sqrt((stop - t) / (stop - (t - 1)) / stop) * rnorm(n=n, mean = 0, sd = 1)
  }
  return(data)
}
# --------------------------------------------- #

n <- 10000
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

# bridge parameters
m <- 10
l <- n / m
q <- m / n
p <- 1 / l

# bridge brownian
W <- matrix(0, nrow = n, ncol = ts)
for (i in 1:l){
  a <- 1 + (i - 1) * m
  b <- i * m
  U <- runif(m, 0, 1)
  V <- (i - 1 + U) / l
  W[a:b,ts] <- sqrt(mat) * qnorm(V)
}

W <- bridge(data = W, stop = ts, n = n)

z <- rnorm(ts * n, 0, 1)
z <- matrix(z, nrow = n, ncol = ts)
Z <- matrix(z, nrow = n, ncol = ts) # Importance z
S <- matrix(s, nrow = n, ncol = ts)
B <- matrix(s, nrow = n, ncol = ts) # Bridge stock price
I <- matrix(log(s), nrow = n, ncol = ts) # Importance stock price

prices_cmc <- sds_cmc <- list()
prices_cv  <- sds_cv  <- list()
prices_ss  <- sds_ss  <- list()
prices_is  <- sds_is  <- list()

price_cmc <- sd_cmc <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_cv  <- sd_cv  <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_ss  <- sd_ss  <- matrix(0, nrow = length(strike), ncol = length(sigma))
price_is  <- sd_is  <- matrix(0, nrow = length(strike), ncol = length(sigma))

h <- 1
i <- 1

for (r in rate){
  for (v in sigma){
    for (k in strike){
      
      drift <- (r - 0.5 * v ** 2) * dt
      vol   <- v * sqrt(dt)
      dc    <- exp(-r * mat)
      
      # importance sampling
      y <- 5
      c <- 0
      
      while (abs(y - c) > 0.0001){
        y <- (c + y) / 2
        mu <- rep(vol * (y + k) / y, ts)
        pi <- log(s) + drift + vol * mu
        for (t in 2:ts){
          mu[t] <- mu[t-1] - vol * exp(pi[t-1]) / (ts * y)
          pi[t] <- pi[t-1] + drift + vol * mu[t]
        }
        pi <- mean(exp(pi))
        c  <- max(pi - k, 0)
      }
      
      Z <- sweep(Z, 2, mu, '+')
      
      # Stock prices
      for (t in 2:ts){
        S[,t] <- S[, t-1] * exp(drift + vol * z[,t])
        B[,t] <- s * exp(drift * t + v * W[,t])
        I[,t] <- I[, t-1] + drift + vol * Z[,t]
      }
      
      arith_mean   <- rowMeans(S)
      geom_mean    <- exp(rowMeans(log(S)))
      bridge_mean  <- rowMeans(B)
      bridge_strat <- matrix(bridge_mean, nrow = m)
      is_mean      <- rowMeans(exp(I))
      
      c_a <- dc * pmax(arith_mean - k, 0)
      c_g <- dc * pmax(geom_mean  - k, 0)
      c_s <- dc * pmax(bridge_strat - k, 0)
      c_i <- dc * pmax(is_mean - k, 0) * 
                exp(-rowSums(Z %*% diag(mu)) + 
                      0.5 * c(mu %*% mu))
      
      c_a_bar <- mean(c_a)
      c_g_bar <- c_g_bar_fct(s, k, r, v, mat)
      
      beta <- beta_fct(c_a, c_a_bar, c_g, c_g_bar)
      
      price_cmc[h, i] <- c_a_bar
      sd_cmc[h, i] <- sd(c_a) / sqrt(n)
      
      price_cv[h, i] <- mean(c_a - beta * (c_g - c_g_bar))
      sd_cv[h, i]    <- sd(c_a - beta * (c_g - c_g_bar)) / sqrt(n)
      
      price_ss[h, i] <- sum(p * colMeans(c_s))
      var_ss <- 0
      for (j in 1:l){
        var    <- p ** 2 * var(c_s[,j]) / q
        var_ss <- var_ss + var
      }
      sd_ss[h, i] <- sqrt(var_ss) / sqrt(n)
      
      price_is[h, i] <- mean(c_i)
      sd_is[h, i]    <- sd(c_i) / sqrt(n)
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
  
  prices_ss[[a]] <- price_ss
  rownames(prices_ss[[a]]) <- paste0('k: ', strike)
  colnames(prices_ss[[a]]) <- paste0('v: ', sigma)
  
  sds_ss[[a]] <- sd_ss
  rownames(sds_ss[[a]]) <- paste0('k: ', strike)
  colnames(sds_ss[[a]]) <- paste0('v: ', sigma)
  
  prices_is[[a]] <- price_is
  rownames(prices_is[[a]]) <- paste0('k: ', strike)
  colnames(prices_is[[a]]) <- paste0('v: ', sigma)
  
  sds_is[[a]] <- sd_is
  rownames(sds_is[[a]]) <- paste0('k: ', strike)
  colnames(sds_is[[a]]) <- paste0('v: ', sigma)
}

prices_cmc
prices_cv
prices_ss
prices_is