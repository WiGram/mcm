dev.off()             
rm(list = ls())       
cat("\014")           
options(scipen = 999) 
library(stats4)       
library(ggplot2)      
set.seed(12345)

# ============================================= #
# ===== Asian options ========================= #
# ============================================= #

# ----- Functions ----------------------------- #
# Kemma and Vorst
c_g_bar_fct <- function(s, k, r, v, mat){
  d_star <- 0.5 * (r - v ** 2 / 6) * mat
  d1     <- (log(s / k) + 0.5 * (r + v ** 2 / 6) * mat) / (v * sqrt(mat / 3))
  d2     <- d1 - v * sqrt(mat / 3)
  exp(-r * mat) * (exp(d_star) * s * pnorm(d1) - k * pnorm(d2))
}

# CV beta
beta_fct <- function(c_1, c_1_bar, c_2, c_2_bar){
  sum((c_1 - c_1_bar) * (c_2 - c_2_bar)) /
    sum((c_2 - c_2_bar) ** 2)
}

# Brownian bridge
bridge <- function(data, stop, n){
  for (t in 2:(stop-1)){
    a <- data[,(t-1)]
    b <- data[,stop]
    data[,t] <- ((stop - t) * a + b) / (stop - (t-1)) +
      sqrt((stop - t) / (stop - (t - 1)) / stop) * rnorm(n=n, mean = 0, sd = 1)
  }
  return(data)
}

# Optimal shift
findShift <- function(r, v, k, ts, dt){
  drift <- (r - 0.5 * v ** 2) * dt
  vol   <- v * sqrt(dt)
  
  a <- 0
  b <- 20
  
  while (abs(a - b) > 0.0001){
    y <- (a + b) / 2
    mu <- rep(vol * (y + k) / y, ts)
    pi <- exp(log(s) + drift + vol * mu)
    for (t in 2:ts){
      mu[t] <- mu[t-1] - vol * pi[t-1] / (ts * y)
      pi[t] <- pi[t-1] * exp(drift + vol * mu[t])
    }
    c <- mean(pi) - k - y
    ifelse(c > 0, a <- y, b <- y)
    return(mu)
  }
}

# Stock price
asianStock <- function(r, v, dt, z){
  drift <- (r - 0.5 * v ** 2) * dt
  vol   <- v * sqrt(dt)
  
  stock <- matrix(s, nrow = n, ncol = ts)
  for (t in 2:ts){
    stock[, t] <- stock[, t-1] * exp(drift + vol * z[, t])
  }
  return(stock)
}

emptyMatrix <- function(strike, sigma){
  return(matrix(0, nrow = length(strike), ncol = length(sigma)))
}
# ============================================= #
# ===== Random number generation ============== #
# ============================================= #

# --------------------------------------------- #
asianRandomNumbers <- function(n, ts, seed = 12345){
  set.seed(seed)
  z <- rnorm(ts * n, 0, 1)
  z <- matrix(z, nrow = n, ncol = ts)
  return(z)
}

# --------------------------------------------- #
asianRandomAV <- function(n, ts){
  z <- asianRandomNumbers(n, ts)
  z <- rbind(z, -z)
  return(z)
}

# --------------------------------------------- #
asianRandomSS <- function(n, m, ts, mat){
  l <- n / m
  Z <- matrix(0, nrow = n, ncol = ts * mat)
  for (i in 1:l){
    a <- 1 + (i - 1) * m
    b <- i * m
    U <- runif(m, 0, 1)
    V <- (i - 1 + U) / l
    Z[a:b,ts] <- sqrt(mat) * qnorm(V)
  }
  
  Z <- bridge(Z, stop = ts, n = n)
  return(Z)
}

# --------------------------------------------- #
asianRandomIS <- function(n, mu, ts){
  z  <- asianRandomNumbers(n, ts)
  Z  <- sweep(z, 2, mu, '+')
  return(Z)
}

# --------------------------------------------- #
asianISSS <- function(z, vol, drift, k, ts, s){
  mu <- findShift(r, v, k, ts, dt)
  
  u  <- mu / sqrt(c(mu %*% mu))
  c  <- diag(ts) - u %*% t(u)
  
  Z  <- matrix(0, nrow = n, ncol = ts)
  for (i in 1:l){
    a <- (i - 1) * m
    for (j in 1:m){
      U <- runif(1,0,1)
      V <- (i - 1 + U) / l
      Z[a + j,] <- u * qnorm(V) + 
        c(c %*% rnorm(ts))
    }
  }
  return(Z)
}
# --------------------------------------------- #

# ============================================= #
# ===== Pricing functions ===================== #
# ============================================= #

# --------------------------------------------- #
asianCMC <- function(s, k, r, v, dt, ts, mat, n){
  ptm <- proc.time()
  
  z <- asianRandomNumbers(n, ts)
  
  stock     <- asianStock(r, v, dt, z)
  arithMean <- rowMeans(stock)
  
  c     <- exp(-r * mat) * pmax(arithMean - k, 0)
  
  price <- mean(c)
  se    <- sd(c) / sqrt(n)
  time  <- proc.time() - ptm
  
  return(list(price = price, se = se, time = time))
}

# --------------------------------------------- #
asianAV <- function(s, k, r, v, dt, ts, mat, n){
  ptm <- proc.time()
  
  z <- asianRandomAV(n / 2, ts)
  
  stock     <- asianStock(r, v, dt, z)
  arithMean <- rowMeans(stock)
  
  c <- 0.5 * exp(-r * mat) * (pmax(arithMean[1:(n/2)] - k, 0) + 
                   pmax(arithMean[(1 + n/2):n] - k, 0))
  
  price <- mean(c)
  se    <- sd(c) / sqrt(n)
  time  <- proc.time() - ptm
  
  return(list(price = price, se = se, time = time))
}

# --------------------------------------------- #
asianCV <- function(s, k, r, v, dt, ts, mat, n){
  ptm <- proc.time()
  
  z     <- asianRandomNumbers(n, ts)
  
  stock <- asianStock(r, v, dt, z)
  arithMean <- rowMeans(stock)
  
  geomMean  <- exp(rowMeans(log(stock)))
  geomMu    <- c_g_bar_fct(s, k, r, v, mat)
  
  c_a       <- exp(-r * mat) * pmax(arithMean - k, 0)
  c_g       <- exp(-r * mat) * pmax(geomMean - k,  0)
  c_a_bar   <- mean(c_a)
  
  beta      <- beta_fct(c_a, c_a_bar, c_g, geomMu)
  
  price <- mean(c_a - beta * (c_g - geomMu))
  se    <- sd(c_a - beta * (c_g - geomMu)) / sqrt(n)
  time  <- proc.time() - ptm
  
  return(list(price = price, se = se, time = time))
}

# --------------------------------------------- #
asianIS <- function(s, k, r, v, dt, ts, mat, n){
  ptm <- proc.time()
  
  mu <- findShift(r, v, k, ts, dt)
  
  z <- asianRandomIS(n, mu, ts)
  
  stock     <- asianStock(r, v, dt, z)
  arithMean <- rowMeans(stock)
  
  c <- exp(-r * mat) * pmax(arithMean - k, 0) * 
    exp(-rowSums(z %*% diag(mu)) + 
          0.5 * c(mu %*% mu))
  
  price <- mean(c)
  se    <- sd(c) / sqrt(n)
  time  <- proc.time() - ptm
  
  return(list(price = price, se = se, time = time))
}

# --------------------------------------------- #
asianSS <- function(s, k, r, v, dt, ts, mat, n, m, p){
  l <- n / m
  q <- m / n
  
  ptm <- proc.time()
  
  z <- asianRandomSS(n, m, ts, mat)
  
  drift <- (r - 0.5 * v ** 2) * dt
  vol   <- v * sqrt(dt)
  
  stock <- matrix(s, nrow = n, ncol = ts)
  for (t in 2:ts){
    stock[,t] <- s * exp(drift * t + v * z[,t])
  }
  
  stockMean  <- rowMeans(stock)
  sMeanStrat <- matrix(stockMean, nrow = m)
  c <- exp(-r * mat) * pmax(sMeanStrat - k, 0)
  
  var_ss <- 0
  for (j in 1:l){
    var    <- p ** 2 * var(c[,j]) / q
    var_ss <- var_ss + var
  }
  
  price <- sum(p * colMeans(c))
  se    <- sqrt(var_ss) / sqrt(n)
  time  <- proc.time() - ptm
  
  return(list(price = price, se = se, time = time))
}
# ============================================= #

# ============================================= #
# ===== Parameterisation ====================== #
# ============================================= #

n      <- 10000
s      <- 100
mat    <- 1
ts     <- 500
dt     <- mat / ts
strike <- 110
rate   <- 0.05
sigma  <- 0.30
# strike <- seq(80, 120, 10)
# rate   <- seq(0.03,0.07,0.02)
# sigma  <- seq(0.2,0.40,0.1)

asianCMC(s, 110, 0.05, 0.2, dt, ts, mat, n)
asianAV( s, 110, 0.05, 0.2, dt, ts, mat, n)
asianCV( s, 110, 0.05, 0.2, dt, ts, mat, n)

# bridge parameters
m <- 10
l <- n / m
q <- m / n
p <- m / n

# General matrix initialistion
cmcP <- avP <- cvP <- ssP <- isP <- emptyMatrix(strike, sigma)
cmcS <- avS <- cvS <- ssS <- isS <- emptyMatrix(strike, sigma)
cmcT <- avT <- cvT <- ssT <- isT <- emptyMatrix(strike, sigma)

cmcPs <- avPs <- cvPs <- ssPs <- isPs <- list()
cmcSs <- avSs <- cvSs <- ssSs <- isSs <- list()
cmcTs <- avTs <- cvTs <- ssTs <- isTs <- list()

h <- 1
i <- 1

for (r in rate){
  for (v in sigma){
    for (k in strike){
      
      cmc <- asianCMC(s, k, r, v, dt, ts, mat, n)
      av  <- asianAV( s, k, r, v, dt, ts, mat, n)
      cv  <- asianCV( s, k, r, v, dt, ts, mat, n)
      is  <- asianIS( s, k, r, v, dt, ts, mat, n)
      ss  <- asianSS( s, k, r, v, dt, ts, mat, n, m, p)
      
      cmcP[h, i] <- cmc$price
      cmcS[h, i] <- cmc$se
      cmcT[h, i] <- cmc$time[[3]]
      
      avP[h, i] <- av$price
      avS[h, i] <- av$se
      avT[h, i] <- av$time[[3]]
      
      cvP[h, i] <- cv$price
      cvS[h, i] <- cv$se
      cvT[h, i] <- cv$time[[3]]
      
      isP[h, i] <- is$price
      isS[h, i] <- is$se
      isT[h, i] <- is$time[[3]]
      
      ssP[h, i] <- ss$price
      ssS[h, i] <- ss$se
      ssT[h, i] <- ss$time[[3]]

      h <- h + 1
    }
    h <- 1
    i <- i + 1
  }
  h <- 1
  i <- 1
  
  a <- paste0('r: ', r)
  cmcPs[[a]] <- cmcP
  rownames(cmcPs[[a]]) <- paste0('k: ', strike)
  colnames(cmcPs[[a]]) <- paste0('v: ', sigma)
  
  cmcSs[[a]] <- cmcS
  rownames(cmcSs[[a]]) <- paste0('k: ', strike)
  colnames(cmcSs[[a]]) <- paste0('v: ', sigma)
  
  cmcTs[[a]] <- cmcT
  rownames(cmcTs[[a]]) <- paste0('k: ', strike)
  colnames(cmcTs[[a]]) <- paste0('v: ', sigma)
  
  avPs[[a]] <- avP
  rownames(avPs[[a]]) <- paste0('k: ', strike)
  colnames(avPs[[a]]) <- paste0('v: ', sigma)
  
  avSs[[a]] <- avS
  rownames(avSs[[a]]) <- paste0('k: ', strike)
  colnames(avSs[[a]]) <- paste0('v: ', sigma)
  
  avTs[[a]] <- avT
  rownames(avTs[[a]]) <- paste0('k: ', strike)
  colnames(avTs[[a]]) <- paste0('v: ', sigma)
  
  cvPs[[a]] <- cvP
  rownames(cvPs[[a]]) <- paste0('k: ', strike)
  colnames(cvPs[[a]]) <- paste0('v: ', sigma)
  
  cvSs[[a]] <- cvS
  rownames(cvSs[[a]]) <- paste0('k: ', strike)
  colnames(cvSs[[a]]) <- paste0('v: ', sigma)
  
  cvTs[[a]] <- cvT
  rownames(cvTs[[a]]) <- paste0('k: ', strike)
  colnames(cvTs[[a]]) <- paste0('v: ', sigma)
  
  isPs[[a]] <- isP
  rownames(isPs[[a]]) <- paste0('k: ', strike)
  colnames(isPs[[a]]) <- paste0('v: ', sigma)
  
  isSs[[a]] <- isS
  rownames(isSs[[a]]) <- paste0('k: ', strike)
  colnames(isSs[[a]]) <- paste0('v: ', sigma)
  
  isTs[[a]] <- isT
  rownames(isTs[[a]]) <- paste0('k: ', strike)
  colnames(isTs[[a]]) <- paste0('v: ', sigma)
  
  ssPs[[a]] <- ssP
  rownames(ssPs[[a]]) <- paste0('k: ', strike)
  colnames(ssPs[[a]]) <- paste0('v: ', sigma)
  
  ssSs[[a]] <- ssS
  rownames(ssSs[[a]]) <- paste0('k: ', strike)
  colnames(ssSs[[a]]) <- paste0('v: ', sigma)
  
  ssTs[[a]] <- ssT
  rownames(ssTs[[a]]) <- paste0('k: ', strike)
  colnames(ssTs[[a]]) <- paste0('v: ', sigma)
}

cmcPs$`r: 0.05`
avPs$`r: 0.05`
cvPs$`r: 0.05`
isPs$`r: 0.05`
ssPs$`r: 0.05`