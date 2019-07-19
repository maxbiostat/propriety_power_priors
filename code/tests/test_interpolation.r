source("Gaussian_data.r")
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
    y0 = y_0,
    n0 = N_0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    m0 = mu_0,
    k0 = kappa_0,
    a_0 = x
  )
}
c_a0 <- Vectorize(c_a0)

constant_data <- read.csv("../data/Gaussian_logCA0.csv")

library(mgcv)
fit.gam <- gam(lc_a0 ~ s(a0), data = constant_data)

K <- 1E4
pred_a0s <- seq(0, max(constant_data$a0), length.out = K)

a0_grid <- data.frame(a0 = pred_a0s, lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))

curve(c_a0, 0, 5, xlab = expression(a[0]), ylab = expression(log(c(a[0]))), lwd = 3, lty = 2)
points(a0_grid$a0, a0_grid$lc_pred)

get_approx_lc <- function(x, grid){
  y <- grid[, 1]
  i <- which.min(abs(x-y))
  if(i != 1){
    x1 <- grid[i, 1]
    x2 <- grid[i + 1, 1]
    y1 <- grid[i, 2]
    y2 <- grid[i + 1, 2]
    ans <- y1 + (y2-y1) * (x-x1)/(x2-x1)
  }else{
    ans <- grid[i, 2]
  }
  return(ans)
}
app_f <- function(x) get_approx_lc(x, grid = a0_grid)
app_f <- Vectorize(app_f)

curve(c_a0, 0, 5, xlab = expression(a[0]), ylab = expression(log(c(a[0]))), lwd = 3, lty = 2)
curve(app_f, 0, 5, lwd = 3, col = 2, add = TRUE)
lines(a0_grid$a0, a0_grid$lc_pred, lwd = 2, lty = 2, col = 3)


curve(c_a0, 0, 5, xlab = expression(a[0]),
      ylab = expression(log(c(a[0]))), xlim = c(0, 1), ylim = c(-600, 10) , lwd = 3, lty = 2)
curve(app_f, 0, 5, lwd = 3, col = 2, add = TRUE)
lines(a0_grid$a0, a0_grid$lc_pred, lwd = 2, lty = 2, col = 3)
