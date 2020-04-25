source("power_priors_aux.r")

scenario <- 4
J <- 20

source(paste("analyses_Bernoulli/data_Bernoulli_scenario_", scenario, ".r", sep = ""))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity
unnorm.bern <- stan_model("stan/simple_Bernoulli_posterior_unnormalised.stan")

bb.data.forposterior <- list(
  N0 = N_0,
  y0 = y_0,
  c = cc,
  d = dd,
  nu = nu,
  eta = eta,
  N = N,
  y = y
)

unnorm.posterior.bern <- sampling(unnorm.bern, data = bb.data.forposterior, refresh = 0, iter = 4000)
print(unnorm.posterior.bern, digits_summary = 2)
pairs(unnorm.posterior.bern, pars = c("a_0", "theta"))

### The approximately normalised prior

constant_data <- read.csv(paste("../data/constant_data/Bernoulli_logCA0_scenario_", scenario, "_J=", J, ".csv", sep = ""))
library(mgcv)
fit.gam <- gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)


approx.normalised.model <- stan_model("stan/simple_Bernoulli_posterior_normalised_approximate.stan")

get_posterior_K <- function(K){
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  bb.data.forposterior$pred_grid_x <- a0_grid$a0
  bb.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.bern <- sampling(approx.normalised.model,
                                         data = bb.data.forposterior, refresh = 0, iter = 4000)
  return(approx.norm.posterior.bern)
}

approx.posterior.bern <- get_posterior_K(2e4)

approx.posterior.bern

### Now the normalised prior
normalised.model <- stan_model("stan/simple_Bernoulli_posterior_normalised.stan")

norm.posterior.bern <- sampling(normalised.model,
                                data = bb.data.forposterior, refresh = 0, iter = 4000)
norm.posterior.bern
pairs(norm.posterior.bern, pars = c("a_0", "theta"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.bern, 'a_0')$a_0
a0.approx <- extract(approx.posterior.bern, 'a_0')$a_0
a0.norm <- extract(norm.posterior.bern, 'a_0')$a_0

a0.dt <-  data.frame(a0 = c(a0.unnorm, a0.norm, a0.approx),
                     normalisation = c( rep("none", length(a0.unnorm)),
                                        rep("exact", length(a0.norm)),
                                        rep("approximate", length(a0.approx)) )
)
##
# Eq 8 in Neuenschwander et al. 2009
posterior_a0_Bernoulli <- function(a_0, y0, n0, y, n, cc, dd, eta, nu, log = FALSE){
  term1 <- lgamma(a_0 * n0 + cc + dd) + lgamma(a_0 * y0 + y + cc) + lgamma( a_0 *(n0-y0) + (n - y) + dd)
  term2 <- lgamma(a_0 * y_0 + cc) + lgamma(a_0 * (n0 - y0) + dd ) + lgamma(a_0 * n0 + n + cc + dd)
  term3 <- dbeta(a_0, shape1 = eta, shape2 = nu, log = TRUE)
  ans <- term1 - term2 + term3
  if(!log) ans <- exp(ans)
  return(ans)
}
post_a0 <- function(x) {
  posterior_a0_Bernoulli(a_0 = x, y0 = y_0, n0 = N_0,
                         y = y, n = N, cc = cc, dd = dd,eta = eta, nu = nu)
}
post_a0  <- Vectorize(post_a0) 
Kp <- integrate(post_a0, 0, 1)$value
norm_post_a0 <- function(x) post_a0(x)/Kp
norm_post_a0  <- Vectorize(norm_post_a0)

##

library(ggplot2)

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, eta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  stat_function(fun = norm_post_a0,
                geom = "line", colour = "black", linetype = "solid") + 
  # ggtitle("Bernoulli") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist

ggsave(paste("../figures/a0_posterior_Bernoulli_scenario_", scenario, "_J=", J, ".pdf", sep = ""), a0_dist)

a0_dist_norm <- ggplot(subset(a0.dt,
                              normalisation != "none"),
                       aes(x = normalisation, y = a0, fill = normalisation, colour = normalisation)) +
  geom_boxplot(alpha = .4) +
  ggtitle("Bernoulli") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist_norm

ggsave("../figures/a0_posterior_Bernoulli_normalisation_comparison.pdf", a0_dist_norm)

###

unnorm.theta.dt <- data.frame(theta = extract(unnorm.posterior.bern, 'theta')$theta)
unnorm.theta.dt$normalisation <- "none"

approx.theta.dt <- data.frame(theta = extract(unnorm.posterior.bern, 'theta')$theta)
approx.theta.dt$normalisation <- "approximate"

norm.theta.dt <- data.frame(theta = extract(norm.posterior.bern, 'theta')$theta)
norm.theta.dt$normalisation <- "exact"

par.posteriors <- rbind(unnorm.theta.dt, approx.theta.dt, norm.theta.dt)

a0_star <- 0.05
a_star <- a0_star*bb.data.forposterior$y0 + bb.data.forposterior$c + bb.data.forposterior$y
b_star <- a0_star*(bb.data.forposterior$N0 - bb.data.forposterior$y0) + bb.data.forposterior$d + (bb.data.forposterior$N - bb.data.forposterior$y)

theta_posterior <- ggplot(data = par.posteriors, aes(x = theta, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  stat_function(fun = function(x) dbeta(x, a_star, b_star),
                geom = "line", colour = "black", linetype = "solid") + 
  # geom_vline(xintercept = y/N, linetype = "dashed") +
  scale_x_continuous(expression(theta), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  theme_bw(base_size = 20)

theta_posterior

ggsave(paste("../figures/parameter_posterior_Bernoulli_scenario_", scenario, "_J=", J, ".pdf"), theta_posterior)

theta_posterior.b <- ggplot(data = par.posteriors, aes(x = normalisation,
                                                       y = theta, colour = normalisation, fill = normalisation)) +
  geom_boxplot(alpha = .4) +
  scale_y_continuous(expression(theta), expand = c(0, 0)) +
  theme_bw(base_size = 20)

theta_posterior.b

###

# discrepancy theta
m1t <- mean(approx.theta.dt$theta)
m2t <- mean(norm.theta.dt$theta)
(m1t-m2t)/abs(m2t)

# discrepancy a_0
m1a <- mean(a0.approx)
m2a <- mean(a0.norm)
(m1a-m2a)/abs(m2a)
