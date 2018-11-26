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

# Kemma and Vorst
geomMuFct <- function(s, k, r, v, mat){
  d_star <- 0.5 * (r - v ** 2 / 6) * mat
  d1     <- (log(s / k) + 0.5 * (r + v ** 2 / 6) * mat) / (v * sqrt(mat / 3))
  d2     <- d1 - v * sqrt(mat / 3)
  exp(-r * mat) * (exp(d_star) * s * pnorm(d1) - k * pnorm(d2))
}

# CV beta
betaFct <- function(c, cMean, control, controlMu){
  sum((c - cMean) * (control - controlMu)) /
    sum((control - controlMu) ** 2)
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
asianRandomISSS <- function(n, m, ts, mat, u, c){
  l <- n / m
  Z <- matrix(0, nrow = n, ncol = ts * mat)
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
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
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
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
}

# --------------------------------------------- #
asianCV <- function(s, k, r, v, dt, ts, mat, n){
  ptm <- proc.time()
  
  z     <- asianRandomNumbers(n, ts)
  
  stock     <- asianStock(r, v, dt, z)
  arithMean <- rowMeans(stock)
  
  geomMean  <- exp(rowMeans(log(stock)))
  geomMu    <- geomMuFct(s, k, r, v, mat)
  
  cArith    <- exp(-r * mat) * pmax(arithMean - k, 0)
  cGeom     <- exp(-r * mat) * pmax(geomMean - k,  0)
  
  beta      <- betaFct(cArith, mean(cArith), cGeom, geomMu)
  
  price <- mean(cArith - beta * (cGeom - geomMu))
  se    <- sd(cArith - beta * (cGeom - geomMu)) / sqrt(n)
  time  <- proc.time() - ptm
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
}

# --------------------------------------------- #
asianTwoCV <- function(s, k, r, v, dt, ts, mat, n, CV = 'stock'){
  ptm <- proc.time()
  
  z     <- asianRandomNumbers(n, ts)
  
  stock     <- asianStock(r, v, dt, z)
  sT        <- stock[, ts]
  
  arithMean <- rowMeans(stock)
  cArith    <- exp(-r * mat) * pmax(arithMean - k, 0)
  
  geomMean  <- exp(rowMeans(log(stock)))  
  cGeom     <- exp(-r * mat) * pmax(geomMean - k,  0)
  geomMu    <- geomMuFct(s, k, r, v, mat)
  
  Y  <- cArith
  if (CV == 'stock'){
    stockMu   <- exp(r * mat) * s
    X  <- cbind(1, cGeom, sT)
    
    beta <- solve(t(X) %*% X) %*% (t(X) %*% Y)
    betaGeom <- beta[2] * (cGeom - geomMu)
    betaTwo  <- beta[3] * (sT - stockMu)
  } else {
    euro      <- exp(- r * mat) * pmax(sT - k, 0)
    euroMu    <- euroBS(s, k, r, v, mat)
    
    X  <- cbind(1, cGeom, euro)
    
    beta <- solve(t(X) %*% X) %*% (t(X) %*% Y)
    betaGeom <- beta[2] * (cGeom - geomMu)
    betaTwo  <- beta[3] * (euro - euroMu)
  }
  
  price <- mean(cArith - betaGeom - betaTwo)
  se    <- sd(cArith - betaGeom - betaTwo) / sqrt(n)
  time  <- proc.time() - ptm
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
}

# --------------------------------------------- #
asianThreeCV <- function(s, k, r, v, dt, ts, mat, n){
  ptm <- proc.time()
  
  z     <- asianRandomNumbers(n, ts)
  
  stock     <- asianStock(r, v, dt, z)                # Control 1
  stockMu   <- exp(r * mat) * s                       # Mu 1
  
  arithMean <- rowMeans(stock)
  cArith    <- exp(-r * mat) * pmax(arithMean - k, 0) # To estimate
  
  sT        <- stock[, ts]
  euro      <- exp(- r * mat) * pmax(sT - k, 0)       # control 2
  euroMu    <- euroBS(s, k, r, v, mat)                # Mu 2

  geomMean  <- exp(rowMeans(log(stock)))  
  cGeom     <- exp(-r * mat) * pmax(geomMean - k,  0) # control 3
  geomMu    <- geomMuFct(s, k, r, v, mat)             # Mu 3
  
  Y  <- cArith
  X  <- cbind(1, cGeom, sT, euro)
  
  beta <- solve(t(X) %*% X) %*% (t(X) %*% Y)
  betaGeom  <- beta[2] * (cGeom - geomMu)
  betaStock <- beta[3] * (sT - stockMu)
  betaEuro  <- beta[4] * (euro - euroMu)
  
  price  <- mean(cArith - betaGeom - betaStock - betaEuro)
  se     <- sd(cArith - betaGeom - betaStock - betaEuro) / sqrt(n)
  time   <- proc.time() - ptm
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
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
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
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
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
}

# --------------------------------------------- #
asianISSS <- function(s, k, r, v, dt, ts, mat, n, m, p){
  l <- n / m
  q <- m / n
  
  ptm <- proc.time()
  
  mu <- findShift(r, v, k, ts, dt)
  u  <- mu / sqrt(c(mu %*% mu))
  c  <- diag(ts)- u %*% t(u)
  
  z <- asianRandomISSS(n, m, ts, mat, u, c)
  
  stock <- asianStock(r, v, dt, z)
  
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
  seTime <- se * time[[3]]
  
  return(list(price = price, se = se, time = time[[3]], seTime = seTime))
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
# strike <- 110
# rate   <- 0.05
# sigma  <- 0.30
strike <- seq(90, 110, 10)
rate   <- seq(0.03,0.05,0.02)
sigma  <- seq(0.2,0.30,0.1)

# bridge parameters
m <- 10
l <- n / m
q <- m / n
p <- m / n

asian <- list()

for (r in rate){
  a <- paste0('r: ', r)
  for (v in sigma){
    b <- paste0('v: ', v)
    for (k in strike){
      c <- paste0('k: ', k)
      
      asian[["cmc"]][[a]][[b]][[c]] <- asianCMC(s, k, r, v, dt, ts, mat, n)
      # asian[["av"]][[a]][[b]][[c]]  <- asianAV( s, k, r, v, dt, ts, mat, n)
      # asian[["cv"]][[a]][[b]][[c]]  <- asianCV( s, k, r, v, dt, ts, mat, n)
      # asian[['cvS']][[a]][[b]][[c]] <- asianTwoCV( s, k, r, v, dt, ts, mat, n, CV = 'stock')
      # asian[['cvE']][[a]][[b]][[c]] <- asianTwoCV( s, k, r, v, dt, ts, mat, n, CV = 'euro')
      # asian[['cv3']][[a]][[b]][[c]] <- asianThreeCV( s, k, r, v, dt, ts, mat, n)
      # asian[["is"]][[a]][[b]][[c]]  <- asianIS( s, k, r, v, dt, ts, mat, n)
      # asian[["ss"]][[a]][[b]][[c]]  <- asianSS( s, k, r, v, dt, ts, mat, n, m, p)
      # asian[["isss"]][[a]][[b]][[c]]  <- asianISSS( s, k, r, v, dt, ts, mat, n, m, p)
    }
    asian[["cmc"]][[a]][[b]]  <- as.data.frame(do.call(cbind, asian[["cmc"]][[a]][[b]] ))
    asian[["av"]][[a]][[b]]   <- as.data.frame(do.call(cbind, asian[["av"]][[a]][[b]]  ))
    asian[["cv"]][[a]][[b]]   <- as.data.frame(do.call(cbind, asian[["cv"]][[a]][[b]]  ))
    asian[["cvS"]][[a]][[b]]  <- as.data.frame(do.call(cbind, asian[["cvS"]][[a]][[b]] ))
    asian[["cvE"]][[a]][[b]]  <- as.data.frame(do.call(cbind, asian[["cvE"]][[a]][[b]] ))
    asian[["cv3"]][[a]][[b]]  <- as.data.frame(do.call(cbind, asian[["cv3"]][[a]][[b]] ))
    asian[["is"]][[a]][[b]]   <- as.data.frame(do.call(cbind, asian[["is"]][[a]][[b]]  ))
    asian[["ss"]][[a]][[b]]   <- as.data.frame(do.call(cbind, asian[["ss"]][[a]][[b]]  ))
    asian[["isss"]][[a]][[b]] <- as.data.frame(do.call(cbind, asian[["isss"]][[a]][[b]]))
  }
}
# 
# cmcPs$`r: 0.05`
# avPs$`r: 0.05`
# cvPs$`r: 0.05`
# isPs$`r: 0.05`
# ssPs$`r: 0.05`