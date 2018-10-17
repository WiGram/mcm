source("C:/Users/wigr11ab/Dropbox/KU/K3/MCM/gbm_function.R")

spot       <- 100
interest   <- 0.01
volatility <- 0.2
mat        <- 10
paths      <- 10
periods    <- 10000

x <- gbm(spot, interest, volatility, mat, paths, periods, anti = TRUE, seed = 0)

plot(x[4,], type = 'l')
