library(ggplot2)
library(reshape)

# ===== Functions ============================= #
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
# ============================================= #

# Brownian Bridge with stratified maturity
n   <- 20
m   <- 2
l   <- n / m
mat <- 1
ts  <- 20
dt  <- mat / ts
W   <- matrix(0, n, ts)

# for each stratum ...
for (i in 1:l){
  a <- 1 + (i - 1) * m
  b <- i * m
  U <- runif(m, 0, 1)
  V <- (i - 1 + U) / l
  W[a:b, ts] <- sqrt(mat) * qnorm(V)
}

W <- bridge(data = W, stop = ts, n = n)

df <- as.data.frame(t(W))
df <- melt(df)
df$rowid <- 1:ts

ggplot(data = df, 
       aes(x = rowid, 
           y = value, 
           group = factor(variable),
           color = factor(variable)),
       guides = FALSE) +
  geom_line() +
  labs(x = '', y = 'Brownian Bridge') +
  theme(legend.position = 'none') + 
  scale_x_continuous(breaks = seq(0,20, 10))

# Comparison with cmc brownian motion
z <- cbind(0, matrix(rnorm(n * (ts-1)), n, ts-1))

for (t in 2:ts){
  z[,t] <- z[,t-1] + sqrt(dt) * z[,t]
}

df_cmc <- as.data.frame(t(z))
df_cmc <- melt(df_cmc)
df_cmc$rowid <- 1:ts

ggplot(data = df_cmc, 
       aes(x = rowid, 
           y = value, 
           group = factor(variable),
           color = factor(variable))) +
  geom_line() +
  labs(x = '', y = 'Regular Brownian motion') +
  theme(legend.position = 'none') + 
  scale_x_continuous(breaks = seq(0,20, 10))
