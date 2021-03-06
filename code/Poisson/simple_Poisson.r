true.lambda <- 2
N0 <- 200
y0 <- rpois(N0, lambda = true.lambda)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

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
a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density(alpha = .4) +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Poisson") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)
a0_dist
ggsave("../figures/a0_posterior_Poisson.pdf", a0_dist)
###

unnorm.lambda.dt <- data.frame(lambda = extract(unnorm.posterior.poisson, 'lambda')$lambda)
unnorm.lambda.dt$normalisation <- "unnormalised"
norm.lambda.dt <- data.frame(lambda = extract(norm.posterior.poisson, 'lambda')$lambda)
norm.lambda.dt$normalisation <- "normalised"
par.posteriors <- rbind(unnorm.lambda.dt, norm.lambda.dt)

lambda_posterior <- ggplot(data = par.posteriors, aes(x = lambda, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  scale_x_continuous(expression(lambda), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  geom_vline(xintercept = true.lambda, linetype = "dashed") +
  theme_bw(base_size = 20)

lambda_posterior

ggsave("../figures/parameter_posterior_Poisson.pdf", lambda_posterior)

#####################

get_estimates_a0 <- function(x){
  po.data <- list(
    N0 = N0,
    y0 = y0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    a_0 = x
  )
  res <- sampling(compiled.model.prior, data = po.data, refresh = 0)
  return(res)
}
K <- 25
maxA <- 5
a0s <- seq(0.01, maxA, length.out = K)
# a0s <- pracma::logseq(0.01, maxA, n = K)
# a0s <- c(seq(0.05, 1, length.out = K-3), seq(1.2, maxA, length.out = 3))
system.time(
  results <- lapply(a0s, get_estimates_a0)  
)

### Sensitivity to a_0: parameter estimates
mean_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r)$summary[1:2, "mean"])        
  )
) 
mean_pars$a0 <- a0s
round(mean_pars, 3)
#
sd_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r)$summary[1:2, "sd"])        
  )
) 
sd_pars$a0 <- a0s
round(sd_pars, 4)


library(ggplot2)

ggplot(data = mean_pars, aes(x = a0, y = lambda)) +
  geom_line() +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(Mean(lambda)), expand = c(0, 0), limits = c(1.5, 2.2)) +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = true.lambda, linetype = "dashed") +
  NULL

ggplot(data = sd_pars, aes(x = a0, y = lambda)) +
  geom_line()+
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(SD(lambda)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL
