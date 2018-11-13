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
# 
# drift <- (r - 0.5 * v ** 2) * mat
# vol   <- v * sqrt(mat)
# dc    <- exp(-r*mat)
# 
# z <- rnorm(n)
# 
# stock  <- s * exp(drift + vol * z)
# payoff <- pmax(stock - k, 0)

# ggplot(NULL) + 
#   geom_point(aes(x = stock, y = payoff),
#              alpha = 0.5, shape = 1) + 
#   labs(x = 'Option payoff', y = 'Stock price')

# Simulation accuracy plot
n <- 10000
price <- bs(price = s, vol = v, mat = mat, rate = r, strike = k)
price <- rep(price, n)

z <- rnorm(n, 0, 1)

dc    <- exp(-r * mat)
drift <- (r - 0.5 * v ** 2) * mat
vol   <- v * sqrt(mat)

# analytic stock price
stock_price <- exp(r * mat) * s
stock_price <- rep(stock_price, n)

s_price    <- s_mat     <- rep(0, n)

cmc_payoff <- cmc_price <- cmc_se <- rep(0, n)
cv_payoff  <- cv_price  <- cv_se  <- rep(0, n)

# experiment
s_error <- price_adj <- rep(0,n)

for (i in 1:(n)){
  s_mat[i]   <- s * exp(drift + vol * z[i])
  s_price[i] <- sum(s_mat) / i
  
  cmc_payoff[i] <- dc * pmax(s_mat[i] - k, 0)
  cmc_price[i]  <- sum(cmc_payoff) / i
  cmc_se[i]     <- sd(cmc_payoff[1:i]) / sqrt(i)
  
  s_error[i]    <- dc * (stock_price[i] - s_mat[i])
  price_adj[i]  <- sum(s_error) / i
  
  if (i > 1){
    beta <- sum((cmc_payoff[1:i] - cmc_price[i]) * (s_mat[1:i] - s_price[i])) /
      sum((s_mat[1:i] - s_price[i]) ** 2)
    
    cv_price[i] <- sum(cmc_payoff[1:i] - beta * (s_mat[1:i] - stock_price[i])) / i
    cv_se[i]    <-  sd(cmc_payoff[1:i] - beta * (s_mat[1:i] - stock_price[i])) / sqrt(i)
  } else {
    cv_price[i] <- cmc_price[i]
  }
  if (i == 2){cv_se[i-1] <- cv_se[i]}
}

cmc_upper     <- cmc_price + price_adj
cmc_adjusted  <- cmc_price + 0.5 * price_adj
cmc_adj_upper <- price_adj + 0.5 * cmc_price

cv_top <- cv_price + 1.96 * cv_se
cv_bot <- cv_price - 1.96 * cv_se
cmc_top <- cmc_price + 1.96 * cmc_se
cmc_bot <- cmc_price - 1.96 * cmc_se

x <- seq(1, n, 1)
ggplot(NULL) +
  geom_line(aes(x = x,
                y = cmc_price,
                col = 'CMC price')) +
  geom_line(aes(x = x,
                y = cmc_upper,
                col = 'Error adjusted price')) +
  geom_line(aes(x = x,
                y = price,
                col = 'Analytic price')) +
  geom_ribbon(aes(x = x,
                  ymin = cmc_adjusted,
                  ymax = cmc_price),
              fill = 'red',
              alpha = 0.2) +
  geom_ribbon(aes(x = x,
                  ymin = cmc_price,
                  ymax = cmc_upper),
              fill = 'blue',
              alpha = 0.2) +
  ylim(2,10) +
  labs(x = 'Simulations', y = 'Option price') +
  scale_colour_manual(values = c('CMC price' = 'red', 
                                 'Error adjusted price' = 'blue',
                                 'Analytic price' = 'green'),
                      breaks=c('CMC price',
                               'Error adjusted price',
                               'Analytic price')) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())

ggplot(NULL) +
  geom_line(aes(x = x,
                y = s_price,
                col = 'CMC price')) +
  geom_line(aes(x = x,
                y = stock_price,
                col = 'Analytic price')) +
  geom_ribbon(aes(x = x,
                  ymin = s_price,
                  ymax = stock_price),
              fill = 'blue',
              alpha = 0.2) +
  geom_ribbon(aes(x = x,
                  ymin = stock_price,
                  ymax = s_price),
              fill = 'blue',
              alpha = 0.2) +
  ylim(100,110) +
  labs(x = 'Simulations', y = 'Stock price') +
  scale_colour_manual(values = c('CMC price' = 'red', 
                                 'Analytic price' = 'green'),
                      breaks=c('CMC price',
                               'Analytic price')) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())


ggplot(NULL) +
  geom_line(aes(x = x,
                y = cv_price,
                colour = 'Control variate')) +
  geom_ribbon(aes(x = x,
                  ymin = cv_bot,
                  ymax = cv_top),
              fill = 'red',
              alpha = 0.2) +
  geom_line(aes(x = x,
                y = cmc_price,
                colour = 'Regular variate')) +
  geom_ribbon(aes(x = x,
                  ymin = cmc_bot,
                  ymax = cmc_top),
              fill = 'blue',
              alpha = 0.2) +
  geom_line(aes(x = x,
                y = price,
                colour = 'Black Scholes price')) +
  ylim(4,8) +
  labs(x = 'Simulations', y = 'Option price') +
  scale_colour_manual(values = c('Control variate' = 'red', 
                                 'Black Scholes price' = 'green', 
                                 'Regular variate' = 'blue'),
                      breaks=c('Control variate',
                               'Black Scholes price',
                               'Regular variate')) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())

# Values coordinates colours with labels
# breaks controls ordering