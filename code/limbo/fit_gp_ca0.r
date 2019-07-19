true.mu <- -100
true.sigma <- .1
N0 <- 2000
y0 <- rnorm(N0, mean = true.mu, sd = true.sigma)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)


mu_0 <- -99
kappa0 <- 200
alpha_0 <- 5
beta_0 <- 5

###
get_c_a0_gaussian <- function(y0, n0, alpha0, beta0, m0, k0, a_0){
  nstar <- a_0 * n0
  ybar <- mean(y0)
  s <- mean( (y0-ybar)^2 )
  kappa_n <- k0 + nstar  
  alpha_n <- alpha0 + nstar/2
  beta_n <- beta0 + .5 * (nstar * s +  (k0 * nstar * (ybar - m0)^2 )/kappa_n ) 
  ans <- lgamma(alpha_n)-lgamma(alpha0)
  ans <- ans + alpha0 * log(beta0) - alpha_n * log(beta_n)
  ans <- ans + .5 *( log(k0) - log(kappa_n) )-nstar/2 * log(2*pi)
  return(ans)
}
#
c_a0 <- function(x) {
  get_c_a0_gaussian(
    y0 = y0,
    n0 = N0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    m0 = mu_0,
    k0 = kappa0,
    a_0 = x
  )
}
c_a0 <- Vectorize(c_a0)

################
# true.as <- seq(0, 5, length.out = 100)
# ca0.data <- data.frame(a0 = true.as, lc_a0 = sapply(true.as, c_a0))

ca0.data <- read.csv("Gaussian_logCA0.csv")


###########################################

options(mc.cores = 4)



library(brms)
m0 <- brm(lc_a0 ~ gp(a0), data = ca0.data, control = list(adapt_delta = .99, max_treedepth = 20))

plot(m0)
# K <- 9
# maxA <- max(ca0.data$a0) + 2
# forpred <- seq(0, maxA, length.out = K)

plot(forpred, y_opt_pred, type = "l", xlab = expression(a[0]), ylab = expression(log(c(a[0]))))
points(lc_a0 ~ a0, ca0.data, add = TRUE, pch = 16, cex = 2)
lines(forpred, sapply(forpred, c_a0), lwd = 2, col = 2, lty = 2)
