## L(D|theta)^a_0 pi(theta)

alpha <- 10
beta <- 10
yy <- 20
NN <- 50
a_0 <- 0
###########################

power_post <- function(p, y, N, a0){
  lp <- a0 * ( y* log(p) + (N-y) * log(1-p) ) + dbeta(x = p, shape1 = alpha, shape2 = beta, log = TRUE)
  return(exp(lp))  
}
power_post <- Vectorize(power_post)

curve(power_post(p = x, y = yy, N = NN, a0 = a_0))

log_sq_lik <- function(p, y, N){
  ( y* log(p) + (N-y) * log(1-p) )^2
}
curve(log_sq_lik(p = x, y = yy, N = NN))

expect <- function(p, y, N, a0){
  ans <- power_post(p = p, y = y, N = N, a0 = a0) * log_sq_lik(p = p, y = y, N = N)
  return(ans)
}
expect <- Vectorize(expect)
curve(expect(p = x, y = yy, N = NN, a0 = a_0))
integrate(function(x) expect(p = x, y = yy, N = NN, a0 = a_0), 0, 1)

## 1/X^2
one_over_x_sq <- function(p, y, N){
  exp(
    - 2*(y* log(p) + (N-y) * log(1-p)) 
  )
}
## test that code for 1/L(D|theta)^2 is correct
# pp <- .2
# L <- pp^yy * (1-pp)^{NN-yy}
# 1/L^2
# one_over_x_sq(p = pp, y = yy, N = NN)

expect2 <- function(p, y, N, a0){
  ans <- power_post(p = p, y = y, N = N, a0 = a0) * one_over_x_sq(p = p, y = y, N = N)
  return(ans)
}
expect2 <- Vectorize(expect2)

curve(expect2(p = x, y = yy, N = NN, a0 = a_0))

integrate(function(x) expect2(p = x, y = yy, N = NN, a0 = a_0), 0, 1)


## -log(X)/X

minus_logX_over_X <- function(p, y, N){
  logX <- y* log(p) + (N-y) * log(1-p)
  return(
    -logX/exp(logX)
  )
}

expect3 <- function(p, y, N, a0){
  ans <- power_post(p = p, y = y, N = N, a0 = a0) * minus_logX_over_X(p = p, y = y, N = N)
  return(ans)
}
expect3 <- Vectorize(expect3)

curve(expect3(p = x, y = yy, N = NN, a0 = a_0))

integrate(function(x) expect3(p = x, y = yy, N = NN, a0 = a_0), 0, 1)
