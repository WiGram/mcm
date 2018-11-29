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
  exp(-r*mat) * (exp(d_star) * s * pnorm(d1) - k * pnorm(d2))
}

# CV beta
beta_fct <- function(c_1, c_1_bar, c_2, c_2_bar){
  sum((c_1 - c_1_bar) * (c_2 - c_2_bar)) /
    sum((c_2 - c_2_bar) ** 2)
}
# --------------------------------------------- #

# Arithmetic mean

# We want to first simulate the log of the stock
set.seed(12345)
n <- 10000
s <- 100
k <- 110
r <- 0.05
v <- 0.20

mat <- 1        # e.g. four months
ts  <- 500       # time steps, e.g. 22 per month
dt  <- mat / ts

z <- rnorm(ts * n, 0, 1)
z <- matrix(z, nrow = n, ncol = ts)
S <- matrix(s, nrow = n, ncol = ts)

drift <- (r - 0.5 * v ** 2) * dt
vol   <- v * sqrt(dt)

for (t in 2:ts){
  S[, t] <- S[, t-1] * exp(drift + vol * z[, t])
}

arith_mean <- rowMeans(S)
geom_mean  <- exp(rowMeans(log(S)))

c_a <- exp(-r * mat) * pmax(arith_mean - k, 0)
c_g <- exp(-r * mat) * pmax(geom_mean  - k, 0)
c_e <- exp(-r * mat) * pmax(S[,ts] - k, 0)

c_a_bar <- mean(c_a)
c_g_bar <- c_g_bar_fct(s, k, r, v, mat)
c_e_bar <- mean(c_e)

# Analytical solutions
s1  <- 2 * ts + 1
s2  <- 6 * (ts + 1)
sig <- v * sqrt(s1 / s2)
rho <- 0.5 * (r - 0.5 * v**2)
d1  <- (log(s / k) + (rho + sig ** 2) * mat) / ( sig*sqrt(mat))
d2  <- d1 - sig * sqrt(mat)

analytical1 <- exp(-r * mat) * (s * exp((rho + 0.5 * sig**2) * mat) * pnorm(d1) - k * pnorm(d2))

h   <- mat / ts
mu  <- log(s) + (r - 0.5 * v ** 2) * (mat + h) / 2
sig <- v * sqrt(h) * sqrt((2*ts + 1) * (ts + 1) / (6 * ts))
d1  <- (mu - log(k) + sig ** 2) / sig
d2  <- d1 - sig

analytical2 <- exp(-r * mat) * (exp(mu + 0.5 * sig ** 2) * pnorm(d1) - k * pnorm(d2))

# ============================================= #
# ===== Illustrating correlations ============= #
# ============================================= #

# Here we try to immitate Glasserman
ggplot(NULL) + 
  geom_point(aes(x = S[, ts], y = c_e),
             alpha = 0.3, fill = NA, shape = 'o') +
  labs(x = 'Stock price', y = 'European call option')

ggplot(NULL) +
  geom_point(aes(x = S[,ts], y = c_a),
             alpha = 0.3, fill = NA, shape = 'o') +
  labs(x = 'Stock price', y = 'Arithmetic Mean Asian call option')

ggplot(NULL) +
  geom_point(aes(x = c_e, y = c_a),
             alpha = 0.3, fill = NA, shape = 'o') +
  labs(x = 'European Call option', y = 'Arithmetic Mean Asian call option')

ggplot(NULL) +
  geom_point(aes(x = c_g, y = c_a),
             alpha = 0.3, fill = NA, shape = 'o') +
  labs(x = 'Geometric Mean Asian call option', y = 'Arithmetic Mean Asian call option')

# ============================================= #
# ===== Using geom as control variate for ===== #
# ===== arithmetic asian call option      ===== #
# ============================================= #
# 
# beta <- beta_fct(c_a, c_a_bar, c_g, c_g_bar)
# 
# price <- mean(c_a - beta * (c_g - c_g_bar))
# se    <- sd(c_a - beta * (c_g - c_g_bar)) / sqrt(n)
# 
# c(analytical  = analytical1,
#   analytical  = analytical2,
#   geometric   = mean(c_g),
#   k_and_v     = c_g_bar,
#   geom_se     = sd(c_g) / sqrt(n),
#   arithmetic  = c_a_bar,
#   arith_se    = sd(c_a) / sqrt(n),
#   arith_cv    = price,
#   arith_cv_se = se)