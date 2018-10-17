#####################
## The LSM-function
#####################

lsm <- function(sim_matrix, strike, spot, r){

  periods <- length(sim_matrix[1,])
  paths   <- length(sim_matrix[,1])
  
  nperiods <- seq(0,periods-1,1)
  
  npaths   <- seq(1,paths, 1)
  
  # Renove t0 because it is trivial
  disc     <- exp( - r * nperiods[-1])
  
  # Generate matrix for actual option pricing
  P           <- matrix( 0, nrow = paths, ncol = periods)
  colnames(P) <- paste0("Period ", nperiods, sep = "")
  rownames(P) <- paste0("Path ",   npaths,   sep = "")
  
  # Matrix of simulated paths (here manufactured from article example)
  S           <- sim_matrix
  colnames(S) <- paste0("Period ", nperiods, sep = "")
  rownames(S) <- paste0("Path ",   npaths,   sep = "")
  
  # First Cashflow at T is done outside the loop
  P[, periods] <- pmax(strike - S[,periods], 0)
  
  #Now do LSM within a loop, backwards iteration
  for (h in (periods-1):2){
    # Only in the money (itm) paths carry a potential trade-off
    itm <- strike - S[, h] > 0
    
    # Essentially just P[x, y]. Applying matrix(P) avoids a dimension-error when using apply()
    p <- matrix(P[,(h+1): periods], ncol = length(P[1,(h+1):periods]))
    d <- disc[1: (periods - h)]
    
    # Regressand Y is the npv of future paths
    Y <- apply(p, 1, function(x) sum(x * d) )
    
    # Considering only itm paths increases computational efficiency
    Y <- Y[itm]
    
    # Regression on two basis functions and a constant. Could increase to e.g. fifth power
    s_itm  <- S[itm, h]
    s_itm2 <- s_itm**2
    one    <- rep(1, length(s_itm))
    X <- cbind(one, s_itm, s_itm2)
    
    # solve() for inverse; %*% for matrix multiplication
    # We get b equivalently from running lm(Y ~ X)
    b <- solve(t(X) %*% X) %*% (t(X) %*% Y)
    
    exercise_value <- pmax(strike - s_itm, 0)
    exp_cont_value <- X %*% b
    exercise       <- exp_cont_value < exercise_value
    idx            <- names(exercise[exercise == T,])
    
    P[idx, h] <- exercise_value[idx]
    P[idx, (h+1):periods] <- 0
  }
  
  Z <- apply(P[,2:periods], 1, function(x) sum(x * disc) )
  
  return(print(mean(Z)))
}