library(ggplot2)
library(reshape2)

n <- 100000
s <- 100
k <- 110
r <- 0.05
v <- 0.2
mat <- 1
dt  <- 1 / 1000
ts  <- mat / dt

set.seed(12345)
z <- matrix(rnorm(n * ts, 0, 1), n, ts)
s_asian <- matrix(s, n, ts)

drift <- (r - 0.5 * v ** 2) * dt
vol   <- v * sqrt(dt)

for (t in 2:ts){
  s_asian[,t] <- s_asian[,(t-1)] * exp(drift + vol * z[,t])
}

s_avg <- rowMeans(s_asian)
c_each <- exp(-r * mat) * pmax(s_avg - k, 0)
c <- mean(c_each)
sd <- sd(c_each) / sqrt(n)
fraction_pos_payoff <- sum(c_each > 0) / length(c_each)

c
sd
fraction_pos_payoff
# Plotting
df <- as.data.frame(t(s_asian))
df <- melt(df)
df$rowid <- 1:ts

ggplot(data = df,
       aes(x = rowid * mat / ts,
           y = value,
           group = factor(variable))) +
  geom_line() +
  labs(x = 'Time (years)', y = 'Stock price') +
  scale_x_continuous(breaks = seq(0,mat,5),
                     minor_breaks = seq(0,mat,2.5)) +
  ylim(25,250)

x <- seq(25, 250, 0.01)
ggplot(NULL) +
  geom_histogram(aes(x = s_asian[, length(s_asian[1,])], y = ..density..),
                 fill = 'white', col = 'black',
                 bins = 100) +
  geom_line(aes(x = x, y = dlnorm(x, log(s) + drift * ts, vol * sqrt(ts))),
            col = 'red') +
  labs(x = "Stock price") + 
  coord_flip()
