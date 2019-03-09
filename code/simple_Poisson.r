true.lambda <- 2
N0 <- 200
y0 <- rpois(N0, lambda = true.lambda)

library(rstan)
rstan_options(auto_write = TRUE)

options(mc.cores = 4)
rstan:::rstudio_stanc("stan/simple_poisson_prior.stan")

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_poisson_prior.stan")

alpha_0 <- 2
beta_0 <- 2
po.data <- list(
  N0 = N0,
  y0 = y0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = .2
)

prior.poisson <- sampling(compiled.model.prior, data = po.data)
prior.poisson
pairs(prior.poisson)

log.c_a0 <- bridgesampling::bridge_sampler(prior.poisson)

get_c_a0_poisson<- function(y0, n0, alpha0, beta0, a_0){
  logPprime <- sum(lfactorial(y0))
  S <- sum(y0)
  ans <- -a_0 * logPprime + lgamma(a_0 * S + alpha0) - (a_0 * S + alpha0)* log(a_0 *n0 + beta0) + (alpha0*log(beta0) - lgamma(alpha0))
  return(ans)
}
exact.log.c_a0 <- get_c_a0_poisson(
    y0 = po.data$y0,
    n0 = po.data$N0,
    alpha0 = po.data$alpha0,
    beta0 = po.data$beta0,
    a_0 = po.data$a_0
  )
exact.log.c_a0
log.c_a0


anaughts <- seq(0, 2, length.out = 100)
c_a0 <- function(x) {
  get_c_a0_poisson(
    y0 = po.data$y0,
    n0 = po.data$N0,
    alpha0 = po.data$alpha0,
    beta0 = po.data$beta0,
    a_0 = x
  )
} 
MaLs <- unlist(parallel::mclapply(anaughts, c_a0, mc.cores = 4))
plot(anaughts, MaLs, xlab = expression(a[0]), main = "Poisson",
     ylab = expression(log(c(a[0]))), type = "l")

#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity
N <- 100
y <- rpois(N, lambda = true.lambda)

compiled.model.unnorm.posterior <- stan_model("stan/simple_poisson_unnormalised_posterior.stan")

delta <- 1
nu <- 1
po.data.forposterior <- list(
  N0 = N0,
  y0 = y0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  delta = delta,
  nu = nu,
  N = N,
  y = y
)

unnorm.posterior.poisson <- sampling(compiled.model.unnorm.posterior, data = po.data.forposterior)
unnorm.posterior.poisson
pairs(unnorm.posterior.poisson, pars = c("a_0", "lambda"))

### Now the normalised prior
compiled.model.norm.posterior <- stan_model("stan/simple_poisson_normalised_posterior.stan")

norm.posterior.poisson <- sampling(compiled.model.norm.posterior, data = po.data.forposterior)
norm.posterior.poisson
pairs(norm.posterior.poisson, pars = c("a_0", "lambda"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.poisson, 'a_0')$a_0
a0.norm <- extract(norm.posterior.poisson, 'a_0')$a_0

a0.dt <- data.frame(a0 = c(a0.unnorm, a0.norm),
                    normalisation = c( rep("unnormalised", length(a0.unnorm)),
                                       rep("normalised", length(a0.norm)) )
)

library(ggplot2)
ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Poisson") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw()

