# Importance sampling examples
dev.off()             # delete plots
rm(list = ls())       # delete variables
cat("\014")           # clear console
options(scipen = 999) # disable scientific notation
library(reshape2)
library(ggplot2)
# ============================================= #
# ===== Stratified Asian options ============== #
# ============================================= #

#' In this script, we stratify the maturity point
#' of a standard Brownian motion. We then bridge
#' the motion by interpolating between midpoints.

#' Finally, We Use the sample paths to generate
#' stock prices on which we evaluate the option
#' prices, applying the CLT.

# ============================================= #
# ===== Functions ============================= #
# ============================================= #

# Intro function: bridging the bastard
bridge <- function(data, stop, n){
  for (t in 2:(stop-1)){
    a <- data[,(t-1)]
    b <- data[,stop]
    data[,t] <- ((stop - t) * a + b) / (stop - (t-1)) +
      sqrt((stop - t) / (stop - (t - 1)) / stop) * rnorm(n=n, mean = 0, sd = 1)
  }
  return(data)
}

# ============================================= #

n <- 100000

#set.seed(12345)
s <- 100
k <- 110
r <- 0.05
v <- 0.2
mat <- 1
ts  <- stop <- 500

m <- 2      # Simulations per stratum
l <- n / m  # Strata (as many as poss')
q <- m / n  # Sims per stratum relative to all sims
p <- 1 / l  # Prob of being in any given stratum

W <- matrix(0, nrow = n, ncol = ts)

# for each stratum ...
for (i in 1:l){
  a <- 1 + (i - 1) * m
  b <- i * m
  U <- runif(m, 0, 1)
  V <- (i - 1 + U) / l
  W[a:b,ts] <- sqrt(mat) * qnorm(V)
}

W <- bridge(data = W, stop = ts, n = n)

# df <- as.data.frame(t(W))
# df <- melt(df)
# df$rowid <- 1:ts
# 
# ggplot(data = df, 
#        aes(x = rowid, 
#            y = value, 
#            group = factor(variable))) +
#   geom_line() +
#   labs(x = '', y = 'Brownian Bridges')

s_ss <- matrix(s, n, ts)
mu  <- (r - 0.5 * v ** 2) / ts
dc  <- exp(-r * mat)

# test <- s * exp(mu * (2:ts) + v * W[,(2:ts)]) # <- should work!!

for (t in 2:ts){
  s_ss[,t]  <- s * exp(mu * t + v * W[,t])
}

# df_s <- as.data.frame(t(s_ss))
# df_s <- melt(df_s)
# df_s$rowid <- 1:ts
# 
# ggplot(data = df_s, 
#        aes(x = rowid, 
#            y = value, 
#            group = factor(variable))) +
#   geom_line() +
#   labs(x = '', y = 'Stock price')

s_avg <- rowMeans(s_ss)
s_strat <- matrix(s_avg, nrow = m)
c_ss  <- dc * pmax(s_strat - k, 0)
price_ss  <- sum(p * colMeans(c_ss))
price_ss

var_ss <- 0
for (j in 1:l){
  var <- p ** 2 * var(c_ss[,j]) / q
  var_ss <- var_ss + var
}
sd_ss <- sqrt(var_ss) / sqrt(n)
sd_ss
