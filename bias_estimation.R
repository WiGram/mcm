mu <- 10
sigma2 <- 10
n <- 100
m <- 1000

theta_hat_1 <- rep(0, n)
theta_hat_2 <- rep(0, n)

for (i in 1:m){
  x <- rnorm(n, mu, sqrt(sigma2))
  
  theta_hat_1[i] <- mean(x)
  theta_hat_2[i] <- 0.5 * (mean(x) + 5 / n)
}

bias_1 <- mean(theta_hat_1 - mu)
var_1  <- var(theta_hat_1)

# The bias can be larger, if the variance is much
# smaller!
bias_2 <- mean(theta_hat_2 - mu)
var_2  <- var(theta_hat_2)

# Provide the histogram
xfit     <- seq(min(theta_hat_1), 
                max(theta_hat_1), 
                length = length(theta_hat_1))
yfit     <- dnorm(xfit,
                  mean = mean(theta_hat_1),
                  sd = sd(theta_hat_1))

ggplot(NULL, aes(theta_hat_1)) + 
  geom_histogram(col = 'white',
                 fill = 'black',
                 bins = 50,
                 aes(y = ..density..)) +
  geom_line(aes(x = xfit, 
                y = yfit, 
                col = 'red')) +
  guides(col = FALSE)

xfit     <- seq(min(theta_hat_2), 
                max(theta_hat_2), 
                length = length(theta_hat_2))
yfit     <- dnorm(xfit,
                  mean = mean(theta_hat_2),
                  sd = sd(theta_hat_2))


ggplot(NULL, aes(theta_hat_2)) + 
  geom_histogram(col = 'white',
                 fill = 'black',
                 bins = 50,
                 aes(y = ..density..)) +
  geom_line(aes(x = xfit, 
                y = yfit, 
                col = 'red')) +
  guides(col = FALSE)


