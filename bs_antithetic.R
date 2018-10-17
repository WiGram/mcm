dev.off()             # delete plots
rm(list = ls())       # delete variables
cat("\014")           # clear console
options(scipen = 999) # disable scientific notation

# ============================================= #
# ===== Simulating european call prices ======= #
# ============================================= #

#' Here I learnt an important lesson: because
#' of Jensen's inequality, it is important when
#' simulating payoffs, to not take the average
#' of the terminal stock price, and the compute
#' the payoff, but instead to compute the mean
#' of the payoffs on each path. Averaging stock
#' price would lead to a downward bias because
#' of the convexity of the call option.

# ============================================= #
# ===== Analytical Black Scholes price ======== #
# ============================================= #

bs <- function(s0=100, v=0.2,mat=10,r=0.05,k=100){
  d1 <- (log(s0 / k) + (r + 0.5 * v ** 2) * mat)/
    (v * sqrt(mat))
  d2 <- d1 - v * sqrt(mat)
  
  #mean defaults to 0, sd defaults to 1
  n1 <- pnorm(d1)
  n2 <- pnorm(d2)
  dc <- exp(-r * mat)
  
  return(s0 * n1 - dc * k * n2)
}

# ============================================= #
# ===== Crude monte carlo simulation (cmc) ==== #
# ============================================= #

sim_european_call<-function(S_0=100,K=100,sig=0.2,
                            mat=10,r=0.05,n=1000){
  z <- rnorm(n)
  
  dc    <- exp(-r * mat)
  drift <- (r - 0.5 * sig ** 2) * mat
  vol   <- sig * sqrt(mat)
  
  #Simulate the stock
  s_T <- S_0 * exp(drift + vol * z)
  
  #Calculate payoffs
  payoffs <- dc * pmax(s_T-K,0)
  
  #Simulate results and bounds
  price <- sum(payoffs) / n
  sd    <- sqrt(sum((payoffs - price) ** 2) / (n - 1))
  se    <- 1.96 * sd / sqrt(n)
  lower <- price - se
  upper <- price + se
  return(c(Price=price,SE=se,Lower=lower,Upper=upper))
}

# ============================================= #
# ===== Antithetic variates (av) simulation === #
# ============================================= #

sim_european_call_av<-function(S_0=100,K=100,sig=0.2,
                               mat=10,r=0.05,n=1000){
  z <- rnorm(n/2)
  
  dc    <- exp(-r * mat)
  drift <- (r - 0.5 * sig ** 2) * mat
  vol   <- sig * sqrt(mat)
  
  #Simulate the payoffs with both processes (sign on vol)
  sim_payoff_1 <- dc * pmax(S_0 * exp(drift + vol * z) - K, 0)
  sim_payoff_2 <- dc * pmax(S_0 * exp(drift - vol * z) - K, 0)
  sim_payoff   <- 0.5 * (sim_payoff_1 + sim_payoff_2)
  
  #Calculate results and bounds
  price <- sum(sim_payoff) / (n / 2)
  sd    <- sqrt(sum((sim_payoff - price) ** 2) / (n/2-1))
  se    <- 1.96 * sd / sqrt(n / 2)
  lower <- price - se
  upper <- price + se
  return(c(Price=price,SE=se,Lower=lower,Upper=upper))
}

# ============================================= #
# ===== AV without averaging - higher SE ====== #
# ============================================= #

sim_european_call_av2<-function(S_0=100,K=100,sig=0.2,
                               mat=10,r=0.05,n=1000){
  z <- rnorm(n/2)
  
  dc    <- exp(-r * mat)
  drift <- (r - 0.5 * sig ** 2) * mat
  vol   <- sig * sqrt(mat)
  
  #Simulate the payoffs with both processes (sign on vol)
  sim_payoff_1 <- dc * pmax(S_0 * exp(drift + vol * z) - K, 0)
  sim_payoff_2 <- dc * pmax(S_0 * exp(drift - vol * z) - K, 0)
  sim_payoff   <- c(sim_payoff_1, sim_payoff_2)
  
  #Calculate results and bounds
  price <- sum(sim_payoff) / n
  sd    <- sqrt(sum((sim_payoff - price) ** 2) / (n-1))
  se    <- 1.96 * sd / sqrt(n)
  lower <- price - se
  upper <- price + se
  return(c(Price=price,SE=se,Lower=lower,Upper=upper))
}

# ============================================= #
# ===== Comparison of methods ================= #
# ============================================= #

# Crude monte carlo 
ptm <- proc.time()
sim_european_call()
proc.time() - ptm

# Antithetic monte carlo
ptm <- proc.time()
sim_european_call_av()
proc.time() - ptm

# Analytical price
bs()

# ============================================= #
# ===== Conclusion ============================ #
# ============================================= #

#' For low amount of simulations Antithetic fares much
#' better. The precision is better and the SE smaller.
#' From a time perspective the two fare equally.

# ===== FIN =================================== #