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

n <- 10000
s <- 100
k <- 110
r <- 0.05
v <- 0.2
mat <- 1
dt  <- 1 / 1000
ts  <- mat / dt

drift <- (r - 0.5 * v ** 2) * mat
vol   <- v * sqrt(mat)
dc    <- exp(-r*mat)

z <- rnorm(n)

stock  <- s * exp(drift + vol * z)
payoff <- pmax(stock - k, 0)

ggplot(NULL) + 
  geom_point(aes(x = stock, y = payoff),
             alpha = 0.5, shape = 1) + 
  labs(x = 'Option payoff', y = 'Stock price')

# Simulation accuracy plot
n <- 10000
price <- bs(price = s, vol = v, mat = mat, rate = r, strike = k)
price <- rep(price, n)
z <- rnorm(n, 0, 1)

price <- bs(price = s, vol = v, mat = mat, rate = r, strike = k)
price <- rep(price, n)

stock_price <- exp(r * mat) * s
stock_price <- rep(stock_price, n)

cv_price <- cmc_price <- s_price <- rep(0, 2 * n)

s_price <- cmc_price <- cv_price <- rep(0, n)
s_mat <- cmc_payoff  <- rep(0, n)
for (i in 1:(n)){
  s_mat[i]   <- s * exp(drift + vol * z[i])
  s_price[i] <- sum(s_mat) / i
  
  cmc_payoff[i] <- dc * pmax(s_mat[i] - k, 0)
  cmc_price[i]  <- sum(cmc_payoff) / i
  
  if (i > 1){
    beta <- sum((cmc_payoff[1:i] - cmc_price[i]) * (s_mat[1:i] - s_price[i])) /
      sum((s_mat[1:i] - s_price[i]) ** 2)
    
    cv_price[i] <- sum(cmc_payoff[1:i] - beta * (s_mat[1:i] - stock_price[i])) / i
  } else {
    cv_price[i] <- cmc_price[i]
  }
}

ggplot(NULL) +
  geom_line(aes(x = seq(1, n, 1),
                y = cmc_price,
                colour = 'Regular variate')) +
  geom_line(aes(x = seq(1, n, 1),
                y = cv_price,
                colour = 'Control variate')) +
  geom_line(aes(x = seq(1, n, 1),
                y = price,
                colour = 'Black Scholes price')) +
#   geom_line(aes(x = seq(1, n, 1),
#                 y = s_price,
#                 colour = 'Numerical stock price average')) + 
#   geom_line(aes(x = seq(1, n, 1),
#                 y = stock_price,
#                 colour = 'Analytic stock price'))
  ylim(4,8) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())
