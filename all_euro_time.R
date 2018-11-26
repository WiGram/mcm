dev.off()                # delete plots
rm(list = ls())          # delete variables
cat("\014")              # clear console
options(scipen = 999)    # disable scientific notation
options(digits.secs = 6) # 6 decimals for Sys.time()

# ============================================= #
# ===== Assisting functions =================== #
# ============================================= #

# --------------------------------------------- #
# Random numbers
euroRandomNumber <- function(n, seed = 12345){
  set.seed(seed)
  return(rnorm(n = n, mean = 0, sd = 1))
}

euroRandomAV <- function(n, ts){
  z <- euroRandomNumber(n)
  z <- c(z, -z)
  return(z)
}

euroRandomSS <- function(m, l){
  z <- matrix(0, m, l)
  for (j in 1:l){
    U <- runif(m, 0, 1)
    V <- (j - 1 + U) / l
    z[, j] <- qnorm(V)
  }
  return(z)
}
# --------------------------------------------- #


euroStock <- function(s, r, v, mat, z){
  return( s * exp( (r - 0.5 * v ** 2) * mat + v * sqrt(mat) * z))
}

euroPayoff <- function(stock, r, mat, k){
  return( exp( - r * mat) * pmax(stock - k, 0))
}

euroPayoffAV <- function(stock, r, mat, k, n){
  c1 <- euroPayoff(stock[1:(n/2)], r, mat, k)
  c2 <- euroPayoff(stock[(n/2 + 1):n], r, mat, k)
  return( 0.5 * (c1 + c2))
}

betaFct <- function(c, c_bar, s, s_bar){
  sum((c - c_bar) * (s - s_bar)) /
    sum((s - s_bar)**2)
}

findShift <- function(r, v, k, ts, dt){
  drift <- (r - 0.5 * v ** 2) * mat
  vol   <- v * sqrt(mat)
  
  a <- 0
  b <- 60
  
  while (abs(a - b) > 0.0000001){
    y <- (a + b) / 2
    mu <- vol * (y + k) / y
    pi <- exp(log(s) + drift + vol * mu)
    c <- mean(pi) - k - y
    ifelse(c > 0, a <- y, b <- y)
    return(mu)
  }
}

sdSS <- function(c, l, p, q){
  var_ss <- 0
  for (j in 1:l){
    var <- p ** 2 * var(c[,j]) / q
    var_ss <- var_ss + var
  }
  return(sqrt(var_ss))
}

# ============================================= #
# ===== Pricing formulas ====================== #
# ============================================= #

euroBS <- function(s, k, r, v, mat){
  dc <- exp(-r * mat)
  mean <- log(s) + (r - 0.5 * v ** 2) * mat
  vol  <- v * sqrt(mat)
  
  d2 <- (mean - log(k)) / vol
  d1 <- d2 + vol
  n2 <- pnorm(d2)
  n1 <- pnorm(d1)
  return(s * n1 - dc * k * n2)
}

euroCMC <- function(s, k, r, v, mat, ts, n){
  ptm <- Sys.time()
  
  z     <- euroRandomNumber(n)
  stock <- euroStock(s, r, v, mat, z)
  c     <- euroPayoff(stock, r, mat, k)
  
  price <- mean(c)
  se    <- sd(c) / sqrt(n)
  time  <- as.numeric(Sys.time() - ptm)
  seTime <- se * time
  
  return(list(price = price, se = se, time = time, seTime = seTime))
}

euroAV <- function(s, k, r, v, mat, ts, n){
  ptm <- Sys.time()
  
  z     <- euroRandomAV(n / 2)
  stock <- euroStock(s, r, v, mat, z)
  
  c <- euroPayoffAV(stock, r, mat, k, n)
  
  price <- mean(c)
  se    <- sd(c) / sqrt(n)
  time  <- as.numeric(Sys.time() - ptm)
  seTime <- se * time
  
  return(list(price = price, se = se, time = time, seTime = seTime))
}

euroCV <- function(s, k, r, v, mat, ts, n){
  ptm <- Sys.time()
  
  z     <- euroRandomNumber(n)
  stock <- euroStock(s, r, v, mat, z)
  
  c     <- euroPayoff(stock, r, mat, k)
  muS   <- exp(r * mat) * s
  
  beta  <- betaFct(c, mean(c), stock, muS)
  
  price <- mean(c - beta * (stock - muS))
  se    <- sd(c - beta * (stock - muS)) / sqrt(n)
  time  <- as.numeric(Sys.time() - ptm)
  seTime <- se * time
  
  return(list(price = price, se = se, time = time, seTime = seTime))
}

euroIS <- function(s, k, r, v, mat, ts, n){
  dt  <- mat / ts
  
  ptm <- Sys.time()
  
  mu  <- findShift(r, v, k, ts, dt)
  z   <- euroRandomNumber(n) + mu
  
  stock <- euroStock(s, r, v, mat, z)
  c     <- exp( - r * mat) * pmax(stock - k, 0) * 
    exp( - z * mu + 0.5 * mu ** 2)
  
  price <- mean(c)
  se    <- sd(c) / sqrt(n)
  time  <- as.numeric(Sys.time() - ptm)
  seTime <- se * time
  
  return(list(price = price, se = se, time = time, seTime = seTime))
}

euroSS <- function(s, k, r, v, mat, ts, n, m){
  l  <- n / m
  p  <- 1 / l
  q  <- m / n
  
  ptm <- Sys.time()
  
  z <- euroRandomSS(m, l)
  
  stock <- euroStock(s, r, v, mat, z)
  c     <- exp(-r * mat) * pmax(stock - k, 0)
  
  price <- sum(p * colMeans(c))
  se    <- sdSS(c, l, p, q) / sqrt(n)
  time  <- as.numeric(Sys.time() - ptm)
  seTime <- se * time
  
  return(list(price = price, se = se, time = time, seTime = seTime))
}

euroISSS <- function(s, k, r, v, mat, ts, n, m){
  l <- n / m
  p <- 1 / l
  q <- m / n
  
  # script beginning
  ptm <- Sys.time()

  mu    <- findShift(r, v, k, ts, dt)
  
  z     <- euroRandomSS(m, l) + mu
  
  stock <- euroStock(s, r, v, mat, z)
  c     <- exp( - r * mat) * pmax(stock - k, 0) * 
    exp( - z * mu + 0.5 * mu ** 2)
  
  price  <- sum(p * colMeans(c))
  se     <- sdSS(c, l, p, q) / sqrt(n)
  time  <- as.numeric(Sys.time() - ptm)
  seTime <- se * time
  
  return(list(price = price, se = se, time = time, seTime = seTime))
}

# ============================================= #
# ===== Black Scholes ========================= #
# ============================================= #

# We want to first simulate the log of the stock
n <- 100000
s <- 100
k <- 110
r <- 0.05
v <- 0.2

mat <- 1
ts  <- 500

m <- 10

# ============================================= #
# ===== Simulating varying amount of paths ==== #
# ============================================= #

n_sims <- 10 ** (3:6)
len    <- length(n_sims)

euros <- list()

for (i in 1:len){
  n     <- n_sims[i]
  a <- paste0('n: ', n)
  euros[["cmc"]][[a]]  <- euroCMC( s, k, r, v, mat, ts, n)
  euros[["av"]][[a]]   <- euroAV(  s, k, r, v, mat, ts, n)
  euros[["cv"]][[a]]   <- euroCV(  s, k, r, v, mat, ts, n)
  euros[["is"]][[a]]   <- euroIS(  s, k, r, v, mat, ts, n)
  euros[["ss"]][[a]]   <- euroSS(  s, k, r, v, mat, ts, n, m)
  euros[["isss"]][[a]] <- euroISSS(s, k, r, v, mat, ts, n, m)
}

euros[["cmc"]]  <- as.data.frame(do.call(rbind, euros[["cmc"]] ))
euros[["av"]]   <- as.data.frame(do.call(rbind, euros[["av"]]  ))
euros[["cv"]]   <- as.data.frame(do.call(rbind, euros[["cv"]]  ))
euros[["is"]]   <- as.data.frame(do.call(rbind, euros[["is"]]  ))
euros[["ss"]]   <- as.data.frame(do.call(rbind, euros[["ss"]]  ))
euros[["isss"]] <- as.data.frame(do.call(rbind, euros[["isss"]]))