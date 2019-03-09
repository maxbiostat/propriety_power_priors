N0 <- 100
y0 <- 10

library(rstan)
rstan_options(auto_write = TRUE)

options(mc.cores = 4)
rstan:::rstudio_stanc("stan/simple_beta_binomial_prior.stan")

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_beta_binomial_prior.stan")

cc <- 1/2
dd <- 1/2
bb.data <- list(
  N0 = N0,
  y0 = y0,
  c = cc,
  d = dd,
  a_0 = .34
)

prior.bb <- sampling(compiled.model.prior, data = bb.data)
prior.bb

pairs(prior.bb)

log.c_a0 <- bridgesampling::bridge_sampler(prior.bb)

get_c_a0_beta_binomial <- function(y0, n0, cc, dd, a_0){
  ans <- a_0 * lchoose(n0, y0) + lbeta(a_0 * y0 + cc, a_0 *(n0 -y0) + dd) - lbeta(cc, dd)
  return(ans)
}

get_c_a0_beta_binomial(
  y0 = bb.data$y0,
  n0 = bb.data$N0,
  cc = bb.data$c,
  dd = bb.data$d,
  a_0 = bb.data$a_0
)
log.c_a0

anaughts <- seq(0, 10, length.out = 100)
c_a0 <- function(x) {
  get_c_a0_beta_binomial(
      y0 = bb.data$y0,
      n0 = bb.data$N0,
      cc = bb.data$c,
      dd = bb.data$d,
      a_0 = x
    )
} 
MaLs <- sapply(anaughts, c_a0)

plot(anaughts, MaLs, xlab = expression(a[0]), main = "Binary meta analysis",
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

compiled.model.unnorm.posterior <- stan_model("stan/simple_beta_binomial_posterior_unnormalised.stan")

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

unnorm.posterior.bb <- sampling(compiled.model.unnorm.posterior, data = bb.data.forposterior)
unnorm.posterior.bb
pairs(unnorm.posterior.bb, pars = c("a_0", "theta"))

### Now the normalised prior

compiled.model.norm.posterior <- stan_model("stan/simple_beta_binomial_posterior_normalised.stan")

norm.posterior.bb <- sampling(compiled.model.norm.posterior, data = bb.data.forposterior)
norm.posterior.bb
pairs(norm.posterior.bb, pars = c("a_0", "theta", "lp__"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.bb, 'a_0')$a_0
a0.norm <- extract(norm.posterior.bb, 'a_0')$a_0

a0.dt <- data.frame(a0 = c(a0.unnorm, a0.norm),
                    normalisation = c( rep("unnormalised", length(a0.unnorm)),
                                       rep("normalised", length(a0.norm)) )
)

library(ggplot2)
ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Beta-binomial") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw()
