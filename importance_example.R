# Importance sampling examples
dev.off()             # delete plots
rm(list = ls())       # delete variables
cat("\014")           # clear console
options(scipen = 999) # disable scientific notation
library(stats4)       # base package; use for mle()
# install.packages("ggplot2") # if not already installed
library(ggplot2)      # Graphing package (sweet)


# ============================================= #
# ===== Black Scholes ========================= #
# ============================================= #

# We want to first simulate the log of the stock
set.seed(12345)
n <- 100000
s <- 100
k <- 110
r <- 0.05
v <- 0.2
t <- 1

# Useful variables
dc   <- exp(-r*t)
mean <- log(s) + (r - 0.5 * v ** 2) * t
vol  <- v * sqrt(t)

# Intermediate calculations
d2 <- (mean - log(k)) / vol
d1 <- d2 + vol
n2 <- pnorm(d2)
n1 <- pnorm(d1)
dc <- exp(-r*t)
c  <- s * n1 - dc * k * n2

# ============================================= #
# ===== CMC and density/histogram plot ======== #
# ============================================= #

s_cmc <- rnorm(n, mean, vol)
c_cmc <- dc * pmax(exp(s_cmc) - k, 0)

c(mu_hat = mean(c_cmc), sd_cmc = sd(c_cmc) / sqrt(n))

# ---------------- Plot ----------------------- #
# x <- seq(25, 250, 0.01)
# ggplot(NULL) +
#   geom_histogram(aes(x = exp(s_cmc), y = ..density..),
#                  fill = 'white', col = 'black') +
#   geom_line(aes(x = x, y = dlnorm(x, mu, sd)), 
#             col = 'red') + 
#   labs(x = "Stock price")
# 
# ggplot(NULL) +
#   geom_point(aes(x = exp(s_cmc), y = c_cmc),
#              alpha = 0.2, fill = NA, col = 'blue') +
#   labs(x = 'Stock price', y = 'Option price')
# 
# ggplot(NULL) +
#   geom_line(aes(x = x, y = dlnorm(x, mu, sd) * 1000),
#             col = 'red') +
#   geom_line(aes(x = x, y = dlnorm(x, mu + log(k/s), sd) * 1000),
#             col = 'purple') +
#   geom_line(aes(x = x, y = dlnorm(x, mu + log(125/s), sd) * 1000),
#             col = 'green') +
#   geom_line(aes(x = x, y = dlnorm(x, mu + log(150/s), sd) * 1000),
#             col = 'magenta') +
#   geom_point(aes(x = exp(s_cmc), y = c_cmc),
#              alpha = 0.1, fill = NA, col = 'blue', shape = 'o') +
#   labs(x = 'Stock price', y = 'Option price')


# ============================================= #
# ===== Simulating varying amount of paths ==== #
# ============================================= #

n_sims <- 10 ** (3:6)
len    <- length(n_sims)

price_cmc <- sd_cmc <- rep(0, len)

for (i in 1:len){
  n     <- n_sims[i]
  s_cmc <- rnorm(n = n, mean = mean, sd = vol)
  c_cmc <- dc * pmax(exp(s_cmc) - k, 0)
  
  price_cmc[i] <- mean(c_cmc)
  sd_cmc[i]    <- sd(c_cmc) / sqrt(n)
}

# ============================================= #
# ===== Importance Sampling =================== #
# ============================================= #

# We shift the mean using k_star
k_star <- 10*9:18
trials <- length(k_star)
mean_k <- log(k_star) + (r - 0.5 * v ** 2) * t

# Matrix and vector initialisation
w <- C <- S <- log_s <- matrix(0, nrow = n, ncol = trials)
k_price_is <- k_sd_is  <- rep(0, trials)

# Important: f distribution must be the k_star == s
index  <- which(k_star == s)

# All computations are done in this loop
for (i in 1:trials){
  log_s[, i] <- rnorm(n = n, mean = mean_k[i], sd = vol)
  S[, i]  <- exp(log_s[, i])
  C[, i]  <- dc * pmax(S[, i] - k, 0)
  
  f       <- dnorm(log_s[, i], mean = mean_k[index], sd = vol)
  g       <- dnorm(log_s[, i], mean = mean_k[i],     sd = vol)
  w[, i]  <- f / g
  
  k_price_is[i] <- mean(C[, i] * w[, i])
  k_sd_is[i]    <- sd(C[, i] * w[, i]) / sqrt(n)
}

opt <- which.min(k_sd_is)

price_is <- sd_is <- rep(0, len)
for (i in 1:len){
  m <- n_sims[i]
  log_s[, i] <- rnorm(n = m, mean = mean_k[opt], sd = vol)
  S[, i] <- exp(log_s[, i])
  C[, i] <- dc * pmax(S[, i] - k, 0)
  
  f      <- dnorm(log_s[, i], mean = mean_k[index], sd = vol)
  g      <- dnorm(log_s[, i], mean = mean_k[opt],   sd = vol)
  
  w[, i] <- f / g
  
  price_is[i] <- mean(C[, i] * w[, i])
  sd_is[i]    <- sd(C[, i] * w[, i]) / sqrt(m)
}

# ============================================= #
# ===== Stratified Sampling =================== #
# ============================================= #

price_ss <- sd_ss <- rep(0,len)
for (i in 1:len){
  nn <- n_sims[i]
  
  l <- nn / 2   # Strata (as many as poss')
  m <- nn / l   # Simulations per stratum
  q <- m  / nn  # Sims per stratum relative to all sims
  p <- 1  / l   # Equiprobable strata
  
  z <- matrix(0, m, l) 
  for (j in 1:l){
    U <- runif(m, 0, 1)
    V <- (j - 1 + U) / l
    z[, j] <- sqrt(t) * qnorm(V)
  }
  
  # Stratified Sampling price
  s_ss  <- mean + vol * z
  c_ss  <- dc * pmax(exp(s_ss) - k, 0)
  
  price_ss[i]  <- sum(p * colMeans(c_ss))
  
  # Computing standard deviation
  var_ss <- 0
  for (j in 1:l){
    var <- p ** 2 * var(c_ss[,j]) / q
    var_ss <- var_ss + var
  }
  sd_ss[i] <- sqrt(var_ss) / sqrt(nn)
}

# ============================================= #
# ===== Antithetic Variates =================== #
# ============================================= #

price_av <- sd_av <- rep(0,len)
for (i in 1:len){
  m  <- n_sims[i] / 2
  s1 <- rnorm(m, mean = mean, sd = vol)
  s2 <- 2 * mean - s1
  c1 <- dc * pmax(exp(s1) - k, 0)
  c2 <- dc * pmax(exp(s2) - k, 0)
  po <- 0.5 * (c1 + c2)
  price_av[i] <- sum(po) / m
  sd_av[i]    <- sd(po) / sqrt(m)
}

# ============================================= #
# ===== Control Variates ====================== #
# ============================================= #

# W = E[X] + beta (Y - E[Y])
# X is the call option, choose Y to be stock

price_cv <- sd_cv <- rep(0, len)
for (i in 1:len){
  n <- n_sims[i]
  s_cv  <- exp(rnorm(n = n, mean = mean, sd = vol))
  s_bar <- s * exp(r * t)
  c_cv  <- dc * pmax(s_cv - k, 0)
  c_bar <- mean(c_cv)
  beta  <- sum((c_cv - c_bar) * (s_cv - mean(s_cv))) /
            sum((s_cv - mean(s_cv))**2)
  price_cv[i] <- mean(c_cv - beta * (s_cv - s_bar))
  sd_cv[i]    <- sd(c_cv - beta * (s_cv - s_bar)) / sqrt(n)
}

# ============================================= #
# ===== Plotting ============================== #
# ============================================= #

# ggplot(NULL) +
#   geom_histogram(aes(S[, 1]), alpha = 0.2, fill = 'blue', bins = 50) + 
#   geom_histogram(aes(S[, 2]), alpha = 0.2, fill = 'red', bins = 50) +
#   geom_histogram(aes(S[, 3]), alpha = 0.2, fill = 'green', bins = 50) + 
#   geom_histogram(aes(S[, 4]), alpha = 0.2, fill = 'yellow', bins = 50)