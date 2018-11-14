library(ggplot2)

# ===== Initial parameters ==================== #
n   <- 10000
s   <- 100
k   <- 160
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

# ===== Initialisation ======================== #
cmc_payoff <- is_payoff <- rep(0, n)
cmc_price  <- is_price  <- rep(0, n)
cmc_se     <- is_se     <- rep(0, n)
w                       <- rep(0, n)
# ============================================= #

# ============================================= #
# ===== importance sampling =================== #
# ============================================= #

k_star <- 90 + seq(0,5*18,5)
trials <- length(k_star)
mean_k <- log(k_star) + drift

# CMC mean identified
index     <- which(k_star == s)
cmc_drift <- mean_k[index]

# Matrix and vector initialisation
k_price_is <- k_sd_is  <- rep(0, trials)

# Choosing optimal drift
for (i in 1:trials){
  S <- mean_k[i] + vol * z
  C <- dc * pmax(exp(S) - k, 0)
  
  f <- dnorm(S, mean = mean_k[index], sd = vol)
  g <- dnorm(S, mean = mean_k[i],     sd = vol)
  w <- f / g
  
  k_price_is[i] <- mean(C * w)
  k_sd_is[i]    <- sd(C * w) / sqrt(n)
}

# Optimal drift specified
opt       <- which.min(k_sd_is)
opt_drift <- mean_k[opt]

# ============================================= #

# ============================================= #
# ===== Showing convergence of price ========== #
# ============================================= #

for (i in 1:n){
  # Crude monte carlo
  cmc_payoff[i] <- dc * pmax(exp(cmc_drift + vol * z[i]) - k, 0)
  cmc_price[i]  <- sum(cmc_payoff) / i
  cmc_se[i]     <- ifelse(i == 1, 0, sd(cmc_payoff[1:i]) / sqrt(i))
  if(i == 2){cmc_se[i-1] <- cmc_se[i]}
  
  # Importance sampling
  s_mat        <- opt_drift + vol * z[i]
  is_payoff[i] <- dc * pmax(exp(s_mat) - k, 0)
  
  f      <- dnorm(s_mat, mean = cmc_drift, sd = vol)
  g      <- dnorm(s_mat, mean = opt_drift, sd = vol)
  w[i]   <- f / g
  
  is_price[i] <- sum(is_payoff * w) / i
  is_se[i]    <- ifelse(i == 1, 0, sd(is_payoff[1:i] * w[1:i]) / sqrt(i))
  if(i == 2){is_se[i-1] <- is_se[i]}
}

is_top <- is_price + 1.96 * is_se
is_bot <- is_price - 1.96 * is_se
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
                y = is_price,
                colour = 'IS price')) +
  geom_ribbon(aes(x = x,
                  ymin = is_bot,
                  ymax = is_top),
              fill = 'red',
              alpha = 0.2) +
  geom_line(aes(x = x,
                y = price,
                colour = 'Black Scholes price')) +
  labs(x = 'Simulations', y = 'Option price') +
  scale_colour_manual(values = c('IS price' = 'red',
                                 'Black Scholes price' = 'green',
                                 'CMC price' = 'blue'),
                      breaks=c('IS price',
                               'Black Scholes price',
                               'CMC price')) +
  theme(legend.position = 'top',
        legend.direction = 'horizontal',
        legend.title = element_blank())

# Values coordinates colours with labels
# breaks controls ordering