library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_inverse_gaussian_prior.stan")

source("inverse_Gaussian_data.r")
summary(y0)

ig.data <- list(
  N0 = N0,
  y0 = y0,
  alpha_l = a.l,
  beta_l = b.l,
  alpha_m = a.m,
  beta_m = b.m,
  a_0 = .32
)

prior.invgauss <- sampling(compiled.model.prior, data = ig.data)
prior.invgauss
pairs(prior.invgauss, pars = c("mu", "lambda"))

log.c_a0 <- bridgesampling::bridge_sampler(prior.invgauss)
log.c_a0$logml

get_c_a0_invgauss <- function(y0, n0, alpha.l, beta.l, alpha.m, beta.m, a_0){
  invgaus_kernel_forint_2 <- function(u, a, b, c, d, s, r){
    (a/(c * u^2) - b/(c* u) + 1)^-d * u^(s-1) * exp(-r * u)
  }
  invgaus_kernel_forint_2 <- Vectorize(invgaus_kernel_forint_2)
  ###########
  logP <- sum(log(y0))
  S <- sum(y0)
  Sprime <- sum(1/y0)
  C <- a_0*Sprime/2 + beta.l
  D <- (a_0*n0)/2 + alpha.l
  lconst <- -0.5 *n0*a_0*log(2*pi) - (3*a_0/2)*logP + lgamma(D) + alpha.m * log(beta.m) +  alpha.l * log(beta.l) - lgamma(alpha.l) - lgamma(alpha.m)
  kf <- function(x, 
                 a = a_0*S/2,
                 b = a_0*n0,
                 c = C,
                 d = D,
                 s = alpha.m,
                 r = beta.m){
    invgaus_kernel_forint_2(u = x, a = a, b = b, c = c, d = d, s = s, r = r)
  }
  Int <- integrate(kf, 0, Inf)
  ###
  ans <- lconst  -D * log(C) +  log(Int$value)
  return(ans)
}
###
exact.log.c_a0 <- get_c_a0_invgauss(
  y0 = ig.data$y0,
  n0 = ig.data$N0,
  alpha.l = ig.data$alpha_l,
  beta.l = ig.data$beta_l,
  alpha.m = ig.data$alpha_m,
  beta.m = ig.data$beta_m,
  a_0 = ig.data$a_0
)
exact.log.c_a0
log.c_a0

get_c_a0_invgauss_brute_force <- function(x){
  ig.data <- list(
    N0 = N0,
    y0 = y0,
    alpha_l = a.l,
    beta_l = b.l,
    alpha_m = a.m,
    beta_m = b.m,
    a_0 = x
  )
  prior.invgauss <- sampling(compiled.model.prior, data = ig.data, refresh = 0)
  res <- list(
    chain =  prior.invgauss,
    logml = bridgesampling::bridge_sampler(prior.invgauss, silent = TRUE)$logml
  )
  return(res)
}
##
c_a0 <- function(x) {
  get_c_a0_invgauss(
    y0 = ig.data$y0,
    n0 = ig.data$N0,
    alpha.l = ig.data$alpha_l,
    beta.l = ig.data$beta_l,
    alpha.m = ig.data$alpha_m,
    beta.m = ig.data$beta_m,
    a_0 = x
  )
} 
c_a0 <- Vectorize(c_a0)
#############
K <- 10
a0s <- seq(.05, 2, length.out = K)
system.time(
  results <- lapply(a0s, get_c_a0_invgauss_brute_force)  
)
mls <- unlist(lapply(results, function(x) x$logml))

## including the (0, 0) point
a0s <- c(0, a0s)
mls <- c(0, mls)


degree <- 2
maxA <- 6
fit0 <- lm(mls ~ poly(a0s, degree))
predAs <- seq(0, maxA, length.out = 100)
maxy <- max(mls)

fit1 <- mgcv::gam(mls ~ s(a0s, k = degree))

# library(brms)
# fit2 <- brm(mls ~ a0s+ I(a0s^2), data = data.frame(mls, a0s))
# plot(marginal_effects(fit2), points = TRUE)

xx <- 5

predict(fit0, newdata = data.frame(a0s = xx))
predict(fit1, newdata = data.frame(a0s = xx))
c_a0(xx)

plot(a0s, mls, xlab = expression(a[0]), ylab = expression(log(c(a[0]))) )
lines(predAs, predict(fit0, newdata = data.frame(a0s = predAs)), col = 2, lwd = 2)
lines(predAs, predict(fit1, newdata = data.frame(a0s = predAs)), col = 3, lwd = 2)
curve(c_a0, 0.001, maxA, lwd = 2, lty = 2, add = TRUE)
legend(x = "topright", legend = c("True", "Polynomial", "Spline"), col = 1:3, lwd = 2, bty = 'n')

plot(a0s, mls, xlab = expression(a[0]), ylab = expression(log(c(a[0]))),
     xlim = c(0, 1), ylim = c(min(mls), max(c_a0(1), 0) ))
lines(predAs, predict(fit0, newdata = data.frame(a0s = predAs)), col = 2, lwd = 2)
lines(predAs, predict(fit1, newdata = data.frame(a0s = predAs)), col = 3, lwd = 2)
curve(c_a0, 0.001, maxA, lwd = 2, lty = 2, add = TRUE)
legend(x = "topright", legend = c("True", "Polynomial", "Spline"), col = 1:3, lwd = 2, bty = 'n')

### Sensitivity to a_0: parameter estimates
mean_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:2, "mean"])        
  )
) 
mean_pars$a0 <- a0s[-1]
round(mean_pars, 3)
#
sd_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:2, "sd"])        
  )
) 
sd_pars$a0 <- a0s[-1]
round(sd_pars, 4)
