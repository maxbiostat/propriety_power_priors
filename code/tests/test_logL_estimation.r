source("../Gaussian_data.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)


#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity

gs.data.forposterior <- list(
  N0 = N_0,
  y0 = y_0,
  mu0 = mu_0,
  kappa0 = kappa_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  delta = delta,
  nu = nu,
  N = N,
  y = y,
  a_0 = .4
)


original <- stan_model("../stan/simple_gaussian_prior.stan")
new <- stan_model("../stan/simple_gaussian_prior_v2.stan")

original.mcmc <- sampling(original, data = gs.data.forposterior)
new.mcmc <- sampling(new, data = gs.data.forposterior)

original.mcmc
new.mcmc
get_deriv <- function(x){
  gs.data.forposterior$a_0 <- x
  new.mcmc <- sampling(new, data = gs.data.forposterior, refresh = 0)
  logL <- extract(new.mcmc, 'logL')$logL
  return(mean(logL))
} 
get_deriv <- Vectorize(get_deriv)

J <- 15
maxA <- 10
a0s <- c(seq(0.05, 1, length.out = J-3), seq(1.2, maxA, length.out = 3))
system.time(
  ca0.prime <- get_deriv(a0s)
)

plot(a0s, ca0.prime)
plot(a0s, ca0.prime, xlim = c(0, 1))
abline(h = 0, lwd = 2, lty = 2)