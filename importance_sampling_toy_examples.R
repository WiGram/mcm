# Importance sampling examples
dev.off()             # delete plots
rm(list = ls())       # delete variables
cat("\014")           # clear console
options(scipen = 999) # disable scientific notation
library(stats4)       # base package; use for mle()
# install.packages("ggplot2") # if not already installed
library(ggplot2)      # Graphing package (sweet)

# ============================================= #
# ===== Toy example =========================== #
# ============================================= #

# We want to simulate some process defined by:
fx <- function(x){ 10 * exp(- 2 * abs(x - 5) ) }

# We draw n numbers from U(min, max)
n    <- 1000000
xmin <- 0
xmax <- 10
x   <- runif(n = 100000, min = 0, max = 10)

# Uniform density pdf(x) = 1 / (b - a)
phi <- function(min, max){1 / (max - min)}

f_x <- fx(x)

# Analytically, mean = 1, variance = 4
c(mean(f_x), var(f_x))


# ============================================= #
# ===== Plotting pdf ========================== #
# ============================================= #
x   <- seq(0, 10, 0.01)
y   <- rep(0.1, length(x))
f_x <- fx(x) * phi(xmin, xmax)

# ggplot(NULL) +
#   geom_line(aes(x = x, y = f_x), col = 'red') +
#   geom_line(aes(x = x, y = y), col = 'blue') +
#   labs(x = '', y = '')

# Plot against a standard normal distribution
y <- 1 / sqrt(2 * pi) * exp(- 0.5 * (x - 5)^2)

# ggplot(NULL) +
#   geom_line(aes(x = x, y = f_x), col = 'red') +
#   geom_line(aes(x = x, y = y), col = 'blue') +
#   labs(x = '', y = '')

# ============================================= #
# ===== Change of measure ===================== #
# ============================================= #

# We will use the normal distribution as below
mu <- 5
sd <- 1

x   <- rnorm(n = n, mean = mu, sd = sd)
phi <- dunif(x, min = xmin, max = xmax)
q_x <- dnorm(x, mean = mu, sd = sd)
w_x <- phi / q_x
f_x <- fx(x)
y_x <- f_x * w_x
c(mean(y_x), var(y_x))

# Conclusion: the mean is preserved, but the
# variance has been reduced 10-fold.

# ============================================= #
# ===== Some other example ==================== #
# ============================================= #

mean <- 0
sd   <- 2
n    <- 100000
x <- rnorm(n = n, mean = mean, sd = sd)

x2 <- x^2
px <- 0.5 * exp(-abs(x))
qx <- 1 / sqrt(2 * pi * sd ** 2) * 
  exp(-0.5 * ((x - mean)/sd) ** 2)

y <- x2 * px / qx

c(mean(y), var(y))