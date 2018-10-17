set.seed(12345)

# ============================================= #
# ===== Parameters ============================ #
# ============================================= #

# Global parameters
n <- 100000 # Total simulations
s <- 100
k <- 110
r <- 0.05
v <- 0.2
t <- 1

# Useful variables
dc   <- exp(-r*t)
mean <- log(s) + (r - 0.5 * v ** 2) * t
vol  <- v * sqrt(t)

# Stratified sampling parameters
l <- 1000   # Strata
m <- n / l  # Simulations per stratum
q <- m / n  # Sims per stratum relative to all sims
p <- 1 / l  # Equiprobable strata

# ============================================= #
# ===== Stratified Sampling =================== #
# ============================================= #

# Stratified sampling
z <- matrix(0, m, l) 
for (i in 1:l){
  U <- runif(m, 0, 1)
  V <- (i - 1 + U) / l
  z[, i] <- sqrt(t) * qnorm(V)
}

# Stratified Sampling price
s_ss  <- mean + vol * z
c_ss  <- dc * pmax(exp(s_ss) - k, 0)
m_ss  <- sum(p * colMeans(c_ss))

# Computing standard deviation
var_ss <- 0
for (i in 1:l){
  var <- p ** 2 * var(c_ss[,i]) / q
  var_ss <- var_ss + var
}
sd_ss <- sqrt(var_ss) / sqrt(n)

# ============================================= #
# ===== Other option pricing computation ====== #
# ============================================= #

# Crude price
s_cmc  <- rnorm(n, mean, vol)
c_cmc  <- dc * pmax(exp(s_cmc) - k, 0)
m_cmc  <- mean(c_cmc)
sd_cmc <- sd(c_cmc) / sqrt(n)

# Analytical price
d2 <- (mean - log(k)) / vol
d1 <- d2 + vol
n2 <- pnorm(d2)
n1 <- pnorm(d1)
dc <- exp(-r*t)
c  <- s * n1 - dc * k * n2

# Final output
c(mu_cmc = mean(c_cmc), sd_cmc = sd_cmc)
c(mu_ss  = mean(c_ss),  sd_ss  = sd_ss)
c

# ============================================= #
#' In conclusion, stratified sampling reduces 
#' variance by an order of magnitude and some.
# ============================================= #

# ============================================= #
# ===== SS for various sizes of stratum ======= #
# ============================================= #

strata <- c(10, 50, 
            100, 500, 
            1000, 5000, 
            10000, 50000)
len <- length(strata)
m_ss <- sd_ss <- rep(0, len)
h <- 1
for (l in strata){
  # l <- 1000
  m <- n / l  # Simulations per stratum
  q <- m / n  # Sims per stratum relative to all sims
  p <- 1 / l  # Equiprobable strata
  
  # Stratified sampling
  z <- matrix(0, m, l) 
  for (i in 1:l){
    U <- runif(m, 0, 1)
    V <- (i - 1 + U) / l
    z[, i] <- sqrt(t) * qnorm(V)
  }
  
  # Stratified Sampling price
  s_ss  <- mean + vol * z
  c_ss  <- dc * pmax(exp(s_ss) - k, 0)
  m_ss[h]  <- sum(p * colMeans(c_ss))
  
  # Computing standard deviation
  var_ss <- 0
  for (i in 1:l){
    var <- p ** 2 * var(c_ss[,i]) / q
    var_ss <- var_ss + var
  }
  sd_ss[h] <- sqrt(var_ss) / sqrt(n)
  h <- h + 1
}
m_ss
sd_ss
sd_cmc / sd_ss
