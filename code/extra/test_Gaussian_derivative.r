source("../power_priors_aux.r")
source("../data_Gaussian.r")

gs.data <- list(
  N0 = N_0,
  y0 = y_0,
  mu0 = mu_0,
  kappa0 = kappa_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = 1
)

###
get_l_a0_gaussian <- function(y0, n0, alpha0, beta0, m0, k0, a_0){
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

############
l_a0 <- function(x) {
  get_l_a0_gaussian(
    y0 = gs.data$y0,
    n0 = gs.data$N0,
    alpha0 = gs.data$alpha0,
    beta0 = gs.data$beta0,
    m0 = gs.data$mu0,
    k0 = gs.data$kappa0,
    a_0 = x
  )
}
l_a0 <- Vectorize(l_a0)
###
l_a0_p <- function(x) numDeriv::grad(l_a0, x)
l_a0_p <- Vectorize(l_a0_p)
l_a0_pp <- function(x) numDeriv::hessian(l_a0, x)
l_a0_pp <- Vectorize(l_a0_pp)

###

find_zero <- function(){
  obj <- function(x) (l_a0_p(x))^2
  res <- optimize(obj, lower = 0, upper = 10)
  return(res$minimum)
}
find_zero()

lhs <- function(a0){
  digamma(alpha_0 + N_0/2 * a0)
}

rhs <- function(a0){
  ybar <- mean(y_0)
  s <-  sum((y_0-ybar)^2)
  kappa_n <- kappa_0 + N_0*a0 
  ##
  delta <- .5 * ( s  + (kappa_0/ kappa_n) * N_0 *(ybar-mu_0)^2)
  t1 <- (delta *(alpha_0 + N_0/2 * a0) )/(N_0 * (delta*a0  + beta_0))
  t2 <- log(delta*a0  + beta_0)
  t3 <- 4/ (a0 *N_0 + kappa_0)
  t4 <- log(2*pi)
  ans <- t1  + t2 + t3 + t4
  return(ans)
}

curve(lhs, 0, 1, lwd = 3, xlab = expression(a[0]))
curve(rhs, 0, 1, add = TRUE, col = 3, lwd = 3)
abline(v = find_zero(), lty = 2, lwd = 2)

find_zero_2 <- function(){
  obj2 <- function(x) (rhs(x) - lhs(x))^2
  res2 <- optimize(obj2, lower = 0, upper = 10)
  return(res2$minimum)
}

curve(l_a0_p)
abline(h = 0, lty = 2, lwd = 2)
abline(v = find_zero_2(), lty = 2, lwd = 2)

lfn <- function(x) l_a0(x + 1)-l_a0(x)
fn <- function(x) exp(lfn(x))

curve(lfn, 0, 10)
curve(fn, 0, 10)
