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

for (t in 2:ts) {
  s_asian[1, t] <- s_asian[1, (t - 1)] * exp(drift + vol * z_p[t])
  s_asian[2, t] <- s_asian[2, (t - 1)] * exp(drift + vol * z_n[t])
}

ggplot(NULL) +
  geom_line(aes(
    x = seq(dt, 1, dt),
    y = s_asian[1, ],
    colour = 'Regular variate'
  )) +
  geom_line(aes(
    x = seq(dt, 1, dt),
    y = s_asian[2, ],
    colour = 'Antithetic variate'
  )) +
  labs(x = 'Time', y = 'Stock price') +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    legend.title = element_blank()
  )

# European call option price
n <- 10000

price <- bs(price = s, vol = v, mat = mat, rate = r, strike = k)
price <- rep(price, n)

z <- rnorm(n, 0, 1)

dc    <- exp(-r * mat)
drift <- (r - 0.5 * v ** 2) * mat
vol   <- v * sqrt(mat)

cmc_payoff <- av_payoff <- rep(0, n)
cmc_price  <- av_price  <- rep(0, n)
cmc_se     <- av_se     <- rep(0, n)

for (i in 1:n){
  payoff_z     <- dc * pmax(s * exp(drift + vol * z[i]) - k, 0)
  payoff_a     <- dc * pmax(s * exp(drift - vol * z[i]) - k, 0)
  
  av_payoff[i] <- 0.5 * (payoff_z + payoff_a)
  av_price[i]  <- sum(av_payoff) / i
  av_se[i]     <- ifelse(i == 1, 0, sd(av_payoff[1:i]) / sqrt(i))
  if(i == 2){ av_se[i-1] <- av_se[i]}

  cmc_payoff[i] <- payoff_z
  cmc_price[i]  <- sum(cmc_payoff) / i
  cmc_se[i]     <- ifelse(i == 1, 0, sd(cmc_payoff[1:i]) / sqrt(i))
  if(i == 2){cmc_se[i-1] <- cmc_se[i]}
}

av_top  <- av_price  + 1.96 * av_se
av_bot  <- av_price  - 1.96 * av_se
cmc_top <- cmc_price + 1.96 * cmc_se
cmc_bot <- cmc_price - 1.96 * cmc_se

x <- seq(1, n, 1)
ggplot(NULL) +
  geom_line(aes(x = x,
                y = cmc_price,
                colour = 'Regular variate')) +
  geom_line(aes(x = x,
                y = av_price,
                colour = 'Antithetic variate')) +
  geom_ribbon(aes(x = x,
                  ymin = av_bot,
                  ymax = av_top),
              fill = 'red',
              alpha = 0.2) +
  geom_line(aes(x = x,
                y = price,
                colour = 'Black Scholes price')) +
  geom_ribbon(aes(x = x,
                  ymin = cmc_bot,
                  ymax = cmc_top),
              fill = 'blue',
              alpha = 0.2) +
  ylim(4,8) +
  labs(x = 'Simulations', y = 'Option price') +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())
