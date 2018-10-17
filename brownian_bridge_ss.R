# Brownian Bridge with stratified maturity
ts <- 10 # time steps
V  <- rep(0, m)
W  <- matrix(0, m, ts)

# for each strata ...
for (i in 1:m){
  U <- runif(1, 0, 1)
  V[i] <- (i - 1 + U) / m
  W[i, ts] <- sqrt(t) * qnorm(V[i])
  for (j in 2:ts){
    z <- rnorm(1,0,1)
    W[i, j] <- (ts - j) / (ts - (j-1)) * W[i, j - 1] +
      (j - (j-1)) / (ts - (j-1)) * W[i, ts] +
      sqrt((ts - j) * (j - (j-1)) / (ts - (j-1))) * z
  }
}

df <- as.data.frame(t(W))
df <- melt(df)
df$rowid <- 1:ts

ggplot(data = df, 
       aes(x = rowid, 
           y = value, 
           group = factor(variable))) +
  geom_line() +
  labs(x = '', y = 'Stock price')
