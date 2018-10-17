# install.packages('ggplot2')
# install.packages('reshape2')
library(ggplot2)
library(reshape2)

n  <- 10
m  <- 500
dt <- 1 / m
s0 <- 100
r  <- 0.05
v  <- 0.2
t  <- 1

set.seed(12345)
x     <- rnorm(n*m, 0, 1)
x     <- matrix(x, n, m)

log_s <- matrix(log(s0), n, m)
mu    <- (r - 0.5 * v ** 2) * dt
vol   <- v * sqrt(dt)

for (i in 2:m){
  log_s[, i] <- log_s[,(i-1)] + mu + vol * x[, i]
}

df <- as.data.frame(t(exp(log_s)))
df <- melt(df)
df$rowid <- 1:m

ggplot(data = df, 
       aes(x = rowid, 
           y = value, 
           group = factor(variable))) +
  geom_line() +
  labs(x = '', y = 'Stock price')
