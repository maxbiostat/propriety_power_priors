N0 <- 100
y0 <- 10

library(rstan)
rstan_options(auto_write = TRUE)

options(mc.cores = 4)
rstan:::rstudio_stanc("stan/simple_bernoulli_prior.stan")

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_bernoulli_prior.stan")

cc <- 1/2
dd <- 1/2
bb.data <- list(
  N0 = N0,
  y0 = y0,
  c = cc,
  d = dd,
  a_0 = .4
)

prior.bern <- sampling(compiled.model.prior, data = bb.data)
prior.bern

pairs(prior.bern)

log.c_a0 <- bridgesampling::bridge_sampler(prior.bern)

get_c_a0_bernoulli <- function(y0, n0, cc, dd, a_0){
  ans <- lbeta(a_0 * y0 + cc, a_0 *(n0 -y0) + dd)
  return(ans)
}

get_c_a0_bernoulli(
  y0 = bb.data$y0,
  n0 = bb.data$N0,
  cc = bb.data$c,
  dd = bb.data$d,
  a_0 = bb.data$a_0
)
log.c_a0

anaughts <- seq(0, 2, length.out = 100)
c_a0 <- function(x) {
  get_c_a0_bernoulli(
      y0 = bb.data$y0,
      n0 = bb.data$N0,
      cc = bb.data$c,
      dd = bb.data$d,
      a_0 = x
    )
} 
MaLs <- sapply(anaughts, c_a0)

plot(anaughts, MaLs, xlab = expression(a[0]), main = "Bernoulli",
     ylab = expression(log(c(a[0]))), type = "l")

plot(anaughts, exp(MaLs), xlab = expression(a[0]),
     ylab = "Marginal likelihood", type = "l")

#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity
N <- 1000
y <- 200

compiled.model.unnorm.posterior <- stan_model("stan/simple_bernoulli_posterior_unnormalised.stan")

delta <- 1
nu <- 1
bb.data.forposterior <- list(
  N0 = N0,
  y0 = y0,
  c = cc,
  d = dd,
  delta = delta,
  nu = nu,
  N = N,
  y = y
)

unnorm.posterior.bern <- sampling(compiled.model.unnorm.posterior, data = bb.data.forposterior)
unnorm.posterior.bern
pairs(unnorm.posterior.bern, pars = c("a_0", "theta"))

### Now the normalised prior

compiled.model.norm.posterior <- stan_model("stan/simple_bernoulli_posterior_normalised.stan")

norm.posterior.bern <- sampling(compiled.model.norm.posterior, data = bb.data.forposterior)
norm.posterior.bern
pairs(norm.posterior.bern, pars = c("a_0", "theta", "lp__"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.bern, 'a_0')$a_0
a0.norm <- extract(norm.posterior.bern, 'a_0')$a_0

a0.dt <- data.frame(a0 = c(a0.unnorm, a0.norm),
                    normalisation = c( rep("unnormalised", length(a0.unnorm)),
                                       rep("normalised", length(a0.norm)) )
)

library(ggplot2)
ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Bernoulli") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw()
