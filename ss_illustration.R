library(ggplot2)

# Initial European call option price function
bs <- function(price, vol, mat, rate, strike) {
  d1 <-
    (log(price / strike) + (rate + 0.5 * vol ** 2) * mat) / (vol * sqrt(mat))
  d2 <- d1 - vol * sqrt(mat)
  
  # mean and sd do not need to be specified: default to N(0,1)
  n1 <- pnorm(d1, mean = 0, sd = 1)
  n2 <- pnorm(d2, mean = 0, sd = 1)
  dc <- exp(-rate * mat)
  
  return(price * n1 - dc * strike * n2)
}

n <- 10000
s <- 100
k <- 110
r <- 0.05
v <- 0.2
mat <- 1
dt  <- 1 / 1000
ts  <- mat / dt

z <- rnorm(n, 0, 1)

m <- 10     # Simulations per stratum
l <- n / m  # Amount of strata
q <- m  / n # Sims per stratum relative to all sims
p <- 1  / l # Equiprobable strata

z_ss <- matrix(0, m, l) 
for (j in 1:l){
  U <- runif(m, 0, 1)
  V <- (j - 1 + U) / l
  z_ss[, j] <- qnorm(V)
}

z_ss <- as.vector(z_ss)

x <- seq(-4+8/n, 4, 8/n)
ggplot(NULL) +
  geom_histogram(aes(x = z, y = ..density..),
                 fill = 'white', col = 'black',
                 bins = 100) +
  geom_line(aes(x = x, y = dnorm(x, 0, 1)),
            col = 'red') +
  labs(x = "Brownian motion")

ggplot(NULL) +
  geom_histogram(aes(x = z_ss, y = ..density..),
                 fill = 'white', col = 'black',
                 bins = 100) +
  geom_line(aes(x = x, y = dnorm(x, 0, 1)),
            col = 'red') +
  labs(x = "Brownian motion")

drift <- (r - 0.5 * v ** 2) * mat
vol   <- v * sqrt(mat)

S  <- s * exp(drift + vol * z)
SS <- s * exp(drift + vol * z_ss)

C  <- pmax(S - k,  0)
CS <- pmax(SS - k, 0)

x_T <- seq(25, 250, 0.01)
ggplot(NULL) +
  geom_histogram(aes(x = S, y = ..density..),
                 fill = 'white', col = 'black',
                 bins = 100) +
  geom_line(aes(x = x_T, y = dlnorm(x_T, log(s) + drift, vol)),
            col = 'red') +
  labs(x = "Terminal stock price (CMC)")

ggplot(NULL) +
  geom_histogram(aes(x = SS, y = ..density..),
                 fill = 'white', col = 'black',
                 bins = 100) +
  geom_line(aes(x = x_T, y = dlnorm(x_T, log(s) + drift, vol)),
            col = 'red') +
  labs(x = "Terminal stock price (SS)")
