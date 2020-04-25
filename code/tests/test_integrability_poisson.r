## L(D|theta)^a_0 pi(theta)

alpha <- 1
beta <- 1
yy <- 200
a_0 <- 0
###########################

power_post <- function(lambda, y, a0){
  lp <- a0 *(y * log(lambda) - lambda * y)  + dgamma(x = lambda, shape = alpha, rate = beta, log = TRUE)
  return(exp(lp))  
}
power_post <- Vectorize(power_post)

curve(power_post(lambda = x, y = yy, a0 = a_0), 0, 10)

log_sq_lik <- function(lambda, y){
  ( y * log(lambda) - lambda * y -lfactorial(y))^2
}
curve(log_sq_lik(lambda = x, y = yy), 0 , 10)

expect <- function(lambda, y,a0){
  ans <- power_post(lambda = lambda, y = y, a0 = a0) * log_sq_lik(lambda = lambda, y = y)
  return(ans)
}
expect <- Vectorize(expect)
curve(expect(lambda = x, y = yy,a0 = a_0), 0, 10)
integrate(function(x) expect(lambda = x, y = yy, a0 = a_0), 0, Inf)

## 1/X^2
one_over_x_sq <- function(lambda, y){
  exp(
    - 2*(y * log(lambda) - lambda * y) 
  )
}
## test that code for 1/L(D|theta)^2 is correct
# pp <- .2
# L <- pp^yy * (1-pp)^{NN-yy}
# 1/L^2
# one_over_x_sq(p = pp, y = yy, N = NN)

expect2 <- function(lambda, y, a0){
  ans <- power_post(lambda = lambda, y = y, a0 = a0) * one_over_x_sq(lambda = lambda, y = y)
  return(ans)
}
expect2 <- Vectorize(expect2)

curve(expect2(lambda = x, y = yy, a0 = a_0), 0, 10)

integrate(function(x) expect2(lambda = x, y = yy, a0 = a_0), 0, Inf)


## -log(X)/X

minus_logX_over_X <- function(lambda, y){
  logX <- y * log(lambda) - lambda * y
  return(
    -logX/exp(logX)
  )
}

expect3 <- function(lambda, y, a0){
  ans <- power_post(lambda = lambda, y = y, a0 = a0) * minus_logX_over_X(lambda = lambda, y = y)
  return(ans)
}
expect3 <- Vectorize(expect3)

curve(expect3(lambda =  x, y = yy, a0 = a_0), 0, 10)

integrate(function(x) expect3(lambda = x, y = yy, a0 = a_0), 0, Inf)
