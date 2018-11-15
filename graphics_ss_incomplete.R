library(ggplot2)

# ===== Initial parameters ==================== #
n   <- 10000
s   <- 100
k   <- 110
r   <- 0.05
v   <- 0.2
mat <- 1
dt  <- 1 / 1000
ts  <- mat / dt
z   <- rnorm(n, 0, 1)

# ===== Useful computations =================== #
dc    <- exp(-r * mat)
drift <- (r - 0.5 * v ** 2) * mat
vol   <- v * sqrt(mat)
# ============================================= #

# ============================================= #
# ===== Analytic European call option price === #
# ============================================= #
bs <- function(price, vol, mat, rate, strike){
  d1 <- (log(price / strike) + (rate + 0.5 * vol ** 2) * mat)/(vol * sqrt(mat))
  d2 <- d1 - vol * sqrt(mat)
  
  # mean and sd do not need to be specified: default to N(0,1)
  n1 <- pnorm(d1, mean = 0, sd = 1)
  n2 <- pnorm(d2, mean = 0, sd = 1)
  dc <- exp(-rate * mat)
  
  return(price * n1 - dc * strike * n2)
}

price <- bs(price = s, vol = v, mat = mat, rate = r, strike = k)
price <- rep(price, n)
# ============================================= #

# ============================================= #
# ===== Stratified Sampling =================== #
# ============================================= #

m   <- 4     # Simulations per stratum
l   <- n / m # Strata
q   <- m / n # = p
p   <- 1 / l # Equiprobable strata

w   <- rep(0, n)
for (i in 1:l){
  a <- 1 + (i - 1) * m
  b <- i * m
  U     <- runif(m, 0, 1)
  V     <- (i - 1 + U) / l
  w[a:b] <- qnorm(V)
}

# Shuffle, so the entries become random
w <- sample(w, length(w))
# ============================================= #

# ============================================= #
# ===== Showing convergence of price ========== #
# ============================================= #

# ===== Initialisation ======================== #
cmc_payoff <- ss_payoff <- rep(0, n)
cmc_price  <- ss_price  <- rep(0, n)
cmc_se     <- ss_se     <- rep(0, n)

# ============================================= #

for (i in 1:n){
  # Crude monte carlo
  cmc_payoff[i] <- dc * pmax(s * exp(drift + vol * z[i]) - k, 0)
  cmc_price[i]  <- sum(cmc_payoff) / i
  cmc_se[i]     <- ifelse(i == 1, 0, sd(cmc_payoff[1:i]) / sqrt(i))
  if(i == 2){cmc_se[i-1] <- cmc_se[i]}
  
  # Stratified sampling
  if(i <= 4){m <- l <- 1}
  ss_payoff[i]  <- dc * pmax(s * exp(drift + vol * w[i]) - k, 0)
  ss_price[i]   <- sum(ss_payoff) / i
  ss_se[i]      <- ifelse(i == 1, 0, sd(ss_payoff[1:i]) / sqrt(i))
  if(i == 2){ss_se[i-1] <- ss_se[i]}
}

ss_top <- ss_price + 1.96 * ss_se
ss_bot <- ss_price - 1.96 * ss_se
cmc_top <- cmc_price + 1.96 * cmc_se
cmc_bot <- cmc_price - 1.96 * cmc_se

x <- seq(1, n, 1)
ggplot(NULL) +
  geom_line(aes(x = x,
                y = cmc_price,
                colour = 'CMC price')) +
  geom_ribbon(aes(x = x,
                  ymin = cmc_bot,
                  ymax = cmc_top),
              fill = 'blue',
              alpha = 0.2) +
  geom_line(aes(x = x,
                y = ss_price,
                colour = 'SS price')) +
  geom_ribbon(aes(x = x,
                  ymin = ss_bot,
                  ymax = ss_top),
              fill = 'red',
              alpha = 0.2) +
  geom_line(aes(x = x,
                y = price,
                colour = 'Black Scholes price')) +
  ylim(4,8) +
  labs(x = 'Simulations', y = 'Option price') +
  scale_colour_manual(values = c('SS price' = 'red', 
                                 'Black Scholes price' = 'green', 
                                 'CMC price' = 'blue'),
                      breaks=c('SS price',
                               'Black Scholes price',
                               'CMC price')) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())

# Values coordinates colours with labels
# breaks controls ordering