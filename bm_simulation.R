# install.packages('ggplot2')
# install.packages('reshape2')
library(ggplot2)
library(reshape2)

n  <- 100
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
  labs(x = '', y = 'Stock price') +
  ylim(25,250)

x <- seq(25, 250, 0.01)
ggplot(NULL) +
  geom_histogram(aes(x = s_asian[, length(s_asian[1,])], y = ..density..),
                 fill = 'white', col = 'black',
                 bins = 100) +
  geom_line(aes(x = x, y = dlnorm(x, log(s) + drift * ts, vol * sqrt(ts))),
            col = 'red') +
  labs(x = "Stock price")
  # + coord_flip()