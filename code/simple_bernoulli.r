source("power_priors_aux.r")
source("Bernoulli_data.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity
unnorm.bern <- stan_model("stan/simple_bernoulli_posterior_unnormalised.stan")

bb.data.forposterior <- list(
  N0 = N_0,
  y0 = y_0,
  c = cc,
  d = dd,
  delta = delta,
  nu = nu,
  N = N,
  y = y
)

unnorm.posterior.bern <- sampling(unnorm.bern, data = bb.data.forposterior, refresh = 0)
unnorm.posterior.bern
pairs(unnorm.posterior.bern, pars = c("a_0", "theta"))

### The approximately normalised prior

constant_data <- read.csv("../data/Bernoulli_logCA0.csv")
library(mgcv)
fit.gam <- gam(lc_a0 ~ s(a0), data = constant_data)

approx.normalised.model <- stan_model("stan/simple_bernoulli_posterior_normalised_approximate.stan")

Ks <- c(50, 100, 1000, 10000)
get_posterior_K <- function(K){
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  bb.data.forposterior$pred_grid_x <- a0_grid$a0
  bb.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.bern <- sampling(approx.normalised.model,
                                         data = bb.data.forposterior, refresh = 0)
  return(approx.norm.posterior.bern)
}

system.time(
  all.approximates <- lapply(Ks, get_posterior_K)  
)

### Now the normalised prior
normalised.model <- stan_model("stan/simple_bernoulli_posterior_normalised.stan")

norm.posterior.bern <- sampling(normalised.model, data = bb.data.forposterior, refresh = 0)
norm.posterior.bern
pairs(norm.posterior.bern, pars = c("a_0", "theta"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.bern, 'a_0')$a_0

a0.approx.list <- lapply(Ks, function(k){
  i <- match(k, Ks)
  a0 <- extract(all.approximates[[i]], 'a_0')$a_0
  return(data.frame(a0 = a0, normalisation = paste("K=", k, sep = "")))
} )

a0.norm <- extract(norm.posterior.bern, 'a_0')$a_0

a0.dt <- 
  rbind(
    data.frame(a0 = c(a0.unnorm, a0.norm),
               normalisation = c( rep("unnormalised", length(a0.unnorm)),
                                  rep("normalised", length(a0.norm)) )
    ),
    do.call(rbind, a0.approx.list)
  )


library(ggplot2)

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Bernoulli") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist

ggsave("../figures/a0_posterior_Bernoulli.pdf", a0_dist)

a0_dist_norm <- ggplot(subset(a0.dt, normalisation != "unnormalised"), aes(x = normalisation, y = a0,
                                                                           fill = normalisation, colour = normalisation)) +
  geom_boxplot(alpha = .4) +
  ggtitle("Bernoulli") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist_norm

ggsave("../figures/a0_posterior_Bernoulli_normalisation_comparison.pdf", a0_dist_norm)

###

unnorm.theta.dt <- data.frame(theta = extract(unnorm.posterior.bern, 'theta')$theta)
unnorm.theta.dt$normalisation <- "unnormalised"

approx.theta.list <- lapply(Ks, function(k) {
  i <- match(k, Ks)
  theta <- extract(all.approximates[[i]], 'theta')$theta
  return(data.frame(theta = theta, normalisation = paste("K=", k, sep = "")))
})

norm.theta.dt <- data.frame(theta = extract(norm.posterior.bern, 'theta')$theta)
norm.theta.dt$normalisation <- "normalised"

par.posteriors <- rbind(unnorm.theta.dt, norm.theta.dt, do.call(rbind, approx.theta.list))

par.posteriors$normalisation <- factor(par.posteriors$normalisation,
                                       levels = c("unnormalised", "normalised", paste("K=", Ks, sep = "")) )

a0_star <- 0.05
a_star <- a0_star*bb.data.forposterior$y0 + bb.data.forposterior$c + bb.data.forposterior$y
b_star <- a0_star*(bb.data.forposterior$N0 - bb.data.forposterior$y0) + bb.data.forposterior$d + (bb.data.forposterior$N - bb.data.forposterior$y)

theta_posterior <- ggplot(data = par.posteriors, aes(x = theta, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  stat_function(fun = function(x) dbeta(x, a_star, b_star),
                geom = "line", colour = "black", linetype = "solid") + 
  scale_x_continuous(expression(theta), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  theme_bw(base_size = 20)

theta_posterior

ggsave("../figures/parameter_posterior_Bernoulli.pdf", theta_posterior)

theta_posterior.b <- ggplot(data = par.posteriors, aes(x = normalisation, y = theta, colour = normalisation, fill = normalisation)) +
  geom_boxplot(alpha = .4) +
  scale_y_continuous(expression(theta), expand = c(0, 0)) +
  theme_bw(base_size = 20)

theta_posterior.b


