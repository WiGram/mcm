# dev.off()             
# rm(list = ls())       
# cat("\014")           
# options(scipen = 999) 
# library(stats4)       
# library(ggplot2)      

n  <- 100000
s  <- 100
k  <- 110
r  <- 0.05
v  <- 0.2
m  <- 1
ts <- 500

asian <- function(n, s, k, r, v, m, ts, mu = 0, seed = 12345){
  set.seed(seed)
  
  dt   <- m / ts
  dc   <- exp(-r*m)
  vol  <- v * sqrt(m)
  z    <- rnorm(n * ts, 0, 1)
  z    <- matrix(z, nrow = n, ncol = ts)
  ifelse(mu == 0,
         z <- z,
         z <- sweep(z, 2, mu, '+')
  )
  S    <- matrix(log(s), nrow = n, ncol = ts)
  
  drift <- (r - 0.5 * v ** 2) * dt
  vol   <- v * sqrt(dt)
  
  for (t in 2:ts){
    S[,t]   <- S[, t-1] + drift + vol * z[,t]
  }
  
  ifelse(mu == 0,
         c <- exp(-r * m) * pmax(rowMeans(exp(S)) - k, 0),
         c <- exp(-r * m) * pmax(rowMeans(exp(S)) - k, 0) * 
           exp(- rowSums(z %*% diag(mu)) + 
                 0.5 * c(mu %*% mu))
  )
  
  c_bar <- mean(c)
  se    <- sd(c) / sqrt(n)
  c(c_bar = c_bar, se = se)
}

find_y <- function(n, s, k, r, v, m, ts, seed = 12345){
  set.seed(seed)
  
  dt    <- m / ts
  drift <- (r - 0.5 * v ** 2) * dt
  vol   <- v * sqrt(dt)
  
  y <- 5
  c <- 0
  
  while (abs(y - c) > 0.0001){
    y <- (c + y) / 2
    Z <- rep(vol * (y + k) / y, ts)
    S <- log(s) + drift + vol * Z
    for (t in 2:ts){
      Z[t] <- Z[t-1]  - vol * exp(S[t-1]) / (ts * y)
      S[t] <- S[t-1] + drift + vol * Z[t]
    }
    S <- mean(exp(S))
    c <- max(S - k, 0)
  }
  
  asian(n, s, k, r, v, m, ts, mu = Z)
  
}

asian(n, s, k, r, v, m, ts)
find_y(n, s, k, r, v, m, ts)

