gbm <- function(spot, interest, volatility, mat, paths, periods, anti = TRUE, seed = 0){
  
  if(seed != 0){
    set.seed(seed)
  }

  time_steps <- periods * mat
  price <- matrix(spot, nrow = paths, ncol = time_steps)
  
  if (anti == TRUE){
    a <- rnorm(n = paths * time_steps / 2)
    z <- matrix( rbind(a, -a), paths, time_steps)
    t <- 1 / (2 * periods)
    
    for (period in 2:time_steps){
      price[,period] <- price[,period - 1] * exp( (interest - 0.5 * volatility**2) * t + 
                                                    volatility * sqrt(t) * z[, period])
    }
  } else {
    a <- rnorm(n = paths * time_steps)
    z <- matrix( a, paths, time_steps)
    t <- 1 / periods
    
    for (period in 2:time_steps){
      price[,period] <- price[,period - 1] * exp( (interest - 0.5 * volatility**2) * t + 
                                                    volatility * sqrt(t) * z[, period])
    }
  }
  
  return(price)
}

gbm(1, 0.1, 0.2, 4, 4, 12345)
