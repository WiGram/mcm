library(ggplot2)

# Initial European call option price function
bs <- function(price, vol, mat, rate, strike){
  d1 <- (log(price / strike) + (rate + 0.5 * vol ** 2) * mat)/(vol * sqrt(mat))
  d2 <- d1 - vol * sqrt(mat)
  
  # mean and sd do not need to be specified: default to N(0,1)
  n1 <- pnorm(d1, mean = 0, sd = 1)
  n2 <- pnorm(d2, mean = 0, sd = 1)
  dc <- exp(-rate * mat)
  
  return(price * n1 - dc * strike * n2)
}

n <- 100
s <- 100
k <- 110
r <- 0.05
v <- 0.2
mat <- 1
dt  <- 1 / 1000
ts  <- mat / dt

z_p <- rnorm(n / dt, 0, 1)
z_n <- -z_p
s_asian <- matrix(s, 2, ts)

drift <- (r - 0.5 * v ** 2) * dt
vol   <- v * sqrt(dt)

for (t in 2:ts){
  s_asian[1,t] <- s_asian[1,(t-1)] * exp(drift + vol * z_p[t])
  s_asian[2,t] <- s_asian[2,(t-1)] * exp(drift + vol * z_n[t])
}

ggplot(NULL) +
  geom_line(aes(x = seq(dt,1,dt),
                y = s_asian[1,],
                colour = 'Regular variate')) +
  geom_line(aes(x = seq(dt,1,dt),
                y = s_asian[2,],
                colour = 'Antithetic variate')) +
  labs(x = 'Time', y = 'Stock price') +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())

# European call option price
n <- 10000

price <- bs(price = s, vol = v, mat = mat, rate = r, strike = k)
price <- rep(price, n)

av_price <- cmc_price <- rep(0, 2 * n)

z <- rnorm(n, 0, 1)

dc    <- exp(-r * mat)
drift <- (r - 0.5 * v ** 2) * mat
vol   <- v * sqrt(mat)

cmc_payoff <- av_payoff <- 0
cmc_price <- av_price <- rep(0, n)
for (i in 1:(n)){
  av_payoff_1  <- dc * pmax(s * exp(drift + vol * z[i]) - k, 0)
  av_payoff_2  <- dc * pmax(s * exp(drift - vol * z[i]) - k, 0)
  av_payoff    <- 0.5 * (av_payoff_1 + av_payoff_2) + av_payoff

  av_price[i]  <- av_payoff / i

  cmc_payoff_1 <- dc * pmax(s * exp(drift + vol * z[i]) - k, 0)
  cmc_payoff   <- cmc_payoff_1 + cmc_payoff

  cmc_price[i] <- cmc_payoff / i
}

ggplot(NULL) +
  geom_line(aes(x = seq(1, n, 1),
                y = cmc_price,
                colour = 'Regular variate')) +
  geom_line(aes(x = seq(1, n, 1),
                y = av_price,
                colour = 'Antithetic variate')) +
  geom_line(aes(x = seq(1, n, 1),
                y = price,
                colour = 'Black Scholes price')) +
  ylim(4,8) +
  labs(x = 'Simulations', y = 'Option price') +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())
