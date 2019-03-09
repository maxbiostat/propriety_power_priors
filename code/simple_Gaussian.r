true.mu <- 1.3
true.sigma <- sqrt(4)
N0 <- 100
y0 <- rnorm(N0, mean = true.mu, sd = true.sigma)

library(rstan)
rstan_options(auto_write = TRUE)

options(mc.cores = 4)
rstan:::rstudio_stanc("stan/simple_gaussian_prior.stan")

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_gaussian_prior.stan")

mu0 <- 0
kappa0 <- 1/2
alpha_0 <- 3
beta_0 <- .5
gs.data <- list(
  N0 = N0,
  y0 = y0,
  mu0 = mu0,
  kappa0 = kappa0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = .94
)

prior.gaussian <- sampling(compiled.model.prior, data = gs.data)
prior.gaussian

pairs(prior.gaussian)

log.c_a0 <- bridgesampling::bridge_sampler(prior.gaussian)

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

get_c_a0_gaussian(
  y0 = gs.data$y0,
  n0 = gs.data$N0,
  alpha0 = gs.data$alpha0,
  beta0 = gs.data$beta0,
  m0 = gs.data$mu0,
  k0 = gs.data$kappa0,
  a_0 = gs.data$a_0
)
log.c_a0

anaughts <- seq(0, 2, length.out = 100)
c_a0 <- function(x) {
  get_c_a0_gaussian(
    y0 = gs.data$y0,
    n0 = gs.data$N0,
    alpha0 = gs.data$alpha0,
    beta0 = gs.data$beta0,
    m0 = gs.data$mu0,
    k0 = gs.data$kappa0,
      a_0 = x
    )
} 
MaLs <- sapply(anaughts, c_a0)

plot(anaughts, MaLs, xlab = expression(a[0]), main = "Gaussian",
     # ylab = expression(log(c(a[0]))), type = "l")

plot(anaughts, exp(MaLs), xlab = expression(a[0]),
     ylab = "Marginal likelihood", type = "l")

#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity
N <- 50
y <- rnorm(N, mean = true.mu, sd = true.sigma)

compiled.model.unnorm.posterior <- stan_model("stan/simple_gaussian_posterior_unnormalised.stan")

delta <- 1
nu <- 1
gs.data.forposterior <- list(
  N0 = N0,
  y0 = y0,
  mu0 = mu0,
  kappa0 = kappa0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  delta = delta,
  nu = nu,
  N = N,
  y = y
)

unnorm.posterior.gaussian <- sampling(compiled.model.unnorm.posterior, data = gs.data.forposterior)
unnorm.posterior.gaussian
pairs(unnorm.posterior.gaussian, pars = c("a_0", "mu", "sigma_sq"))

### Now the normalised prior
compiled.model.norm.posterior <- stan_model("stan/simple_gaussian_posterior_normalised.stan")

norm.posterior.gaussian <- sampling(compiled.model.norm.posterior, data = gs.data.forposterior)
norm.posterior.gaussian
pairs(norm.posterior.gaussian, pars = c("a_0", "mu", "sigma_sq"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.gaussian, 'a_0')$a_0
a0.norm <- extract(norm.posterior.gaussian, 'a_0')$a_0

a0.dt <- data.frame(a0 = c(a0.unnorm, a0.norm),
                    normalisation = c( rep("unnormalised", length(a0.unnorm)),
                                       rep("normalised", length(a0.norm)) )
)

library(ggplot2)
ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Gaussian") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw()
