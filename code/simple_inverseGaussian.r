######## Simple Gaussian distribution example.
######## This script will illustrate the effect of (not) including the normalising constant c(a0).
######## It will also show the effect of various approximations to c(a0)

source("inverse_Gaussian_data.r")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)


#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity

unnormalised.model <- stan_model("stan/simple_inverse_gaussian_posterior_unnormalised.stan")

ig.data.forposterior <- list(
  N0 = N0,
  y0 = y0,
  alpha_l = a.l,
  beta_l = b.l,
  alpha_m = a.m,
  beta_m = b.m,
  delta = delta,
  nu = nu,
  N = N,
  y = y
)

unnorm.posterior.invgaussian <- sampling(unnormalised.model, data = ig.data.forposterior)
unnorm.posterior.invgaussian
pairs(unnorm.posterior.invgaussian, pars = c("a_0", "mu", "lambda"))

### The approximately normalised prior

constant_data <- read.csv("../data/InvGaussian_logCA0.csv")
library(mgcv)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0), data = constant_data)

plot(constant_data$a0, predict.gam(fit.gam), type = "l")
points(lc_a0 ~ a0 , constant_data)

approx.normalised.model <- stan_model("stan/simple_inverse_gaussian_posterior_normalised_approximate.stan")

Ks <- c(50, 100, 1000, 10000)

get_posterior_K <- function(K){
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  ig.data.forposterior$pred_grid_x <- a0_grid$a0
  ig.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.invgauss <- sampling(approx.normalised.model,
                                          data = ig.data.forposterior, refresh = 0)
  return(approx.norm.posterior.invgauss)
}

system.time(
  all.approximates <- lapply(Ks, get_posterior_K)  
)


##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.invgaussian, 'a_0')$a_0

a0.approx.list <- lapply(Ks, function(k){
  i <- match(k, Ks)
  a0 <- extract(all.approximates[[i]], 'a_0')$a_0
  return(data.frame(a0 = a0, normalisation = paste("K=", k, sep = "")))
} )


a0.dt <- 
  rbind(
    data.frame(a0 = a0.unnorm,
               normalisation = rep("unnormalised", length(a0.unnorm))
    ),
    do.call(rbind, a0.approx.list)
  )


library(ggplot2)

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Inverse-Gaussian") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist

ggsave("../figures/a0_posterior_invgaussian.pdf", a0_dist, dpi = 300)

a0_dist_norm <- ggplot(subset(a0.dt, normalisation != "unnormalised"),
                       aes(x = normalisation, y = a0,
                           fill = normalisation, colour = normalisation)) +
  geom_boxplot(alpha = .4) +
  ggtitle("Gaussian") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist_norm

ggsave("../figures/a0_posterior_invgaussian_normalisation_comparison.pdf", a0_dist_norm, dpi = 300)

###
unnorm.pars <- as.data.frame(
  cbind(extract(unnorm.posterior.invgaussian, 'mu')$mu,
        extract(unnorm.posterior.invgaussian, 'lambda')$lambda)
)
names(unnorm.pars) <- c("mu", "lambda")
#

app.norm.pars.list <- lapply(all.approximates, function(x) {
  pars <- as.data.frame(
    cbind(extract(x, 'mu')$mu,
          extract(x, 'lambda')$lambda)
  )
  names(pars) <- c("mu", "lambda")
  return(pars)
})

library(reshape2)
app.norm.pars.df.list <- lapply(seq_along(app.norm.pars.list), function(i) {
  y <- app.norm.pars.list[[i]]
  dt <- melt(y, variable.name = "parameter", value.name = "sample")
  dt$normalisation <- paste("K=", Ks[i] ,sep = "")
  return(dt)
} )


unnorm.dt <- melt(unnorm.pars, variable.name = "parameter", value.name = "sample")
unnorm.dt$normalisation <- "unnormalised"
posterior.dt <- rbind(unnorm.dt, do.call(rbind, app.norm.pars.df.list))

posterior.dt$normalisation <- factor(posterior.dt$normalisation,
                                  levels = c("unnormalised", paste("K=", Ks, sep = "")) )

###

true.pars <- data.frame(parameter = names(unnorm.pars), value = c(true.mu, true.lambda) )

parameter_posteriors <- ggplot(data = posterior.dt,
                               aes(x = sample, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  scale_x_continuous("") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  facet_wrap(.~parameter,  scales = "free") +
  geom_vline(data = true.pars, aes(xintercept = value), linetype = "dashed") +
  theme_bw(base_size = 20)

parameter_posteriors

ggsave(filename = "../figures/parameter_posterior_invgaussian.pdf", plot = parameter_posteriors, dpi = 300)
