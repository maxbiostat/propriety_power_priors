source("data_logistic_regression_multiple_historical.r")
source("../power_priors_aux.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
library(mgcv)
library(reshape2)

############
### Unnormalised posterior

compiled.model.unnorm.posterior <- stan_model("stan/simple_logistic_regression_posterior_unnormalised_multiple_historical.stan")

N0s <- unlist(lapply(dataSets, nrow))
N0max <- max(N0s)
X0s <- array(NA, dim = c(N0max, P, D))
y0s <- matrix(NA, nrow = N0max, ncol = D)
for(d in 1:D){
  y0s[1:N0s[d], d] <- dataSets[[d]][, (P + 1)]
  X0s[1:N0s[d], 1:P, d] <- as.matrix(dataSets[[d]][, -(P + 1)])
}

lgr.data.forposterior <- list(
  D = D,
  N0max = N0max,
  N0 = N0s,
  X0 = X0s,
  P = P,
  y0 = y0s,
  eta = eta,
  nu = nu,
  X = X,
  N = N,
  y = y
)

unnorm.posterior.lgr <- sampling(compiled.model.unnorm.posterior, data = lgr.data.forposterior)
unnorm.posterior.lgr
pairs(unnorm.posterior.lgr, pars = c("a_0", "beta"))

### The approximately normalised prior


constant.data <- lapply(1:D, function(d) 
  read.csv(paste("data/RegressionLogistic_logCA0_J=20_dataset_", d, ".csv", sep = "")) ) 

J <- 20
fit_gam <- function(dt) mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = dt)
fits <- lapply(constant.data, fit_gam)

approx.normalised.model <- stan_model("stan/simple_logistic_regression_posterior_approximate_multiple_historical.stan")

K <- 20000
maxA0 <- 1
pred.a0s <- seq(0, maxA0, length.out = K)

pred.xs <- matrix(NA, nrow = K, ncol = D)
pred.ys <- matrix(NA, nrow = K, ncol = D)
  
for (d in 1:D){
  pred.xs[, d] <- pred.a0s
  pred.ys[, d] <- predict(fits[[d]], newdata = data.frame(a0 = pred.a0s))
}
lgr.data.forposterior$K <- K
lgr.data.forposterior$pred_grid_x <- pred.xs
lgr.data.forposterior$pred_grid_y <- pred.ys
full.approx.time <- system.time(
  approx.norm.posterior.lgr <- sampling(approx.normalised.model,
                                        data = lgr.data.forposterior, refresh = 0)  
)
full.approx.time

approx.norm.posterior.lgr
pairs(approx.norm.posterior.lgr, pars = c("a_0", "beta"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.lgr, 'a_0')$a_0
a0.approx <- extract(approx.norm.posterior.lgr, 'a_0')$a_0

a0s.unnorm.post <- data.frame(a0 = a0.unnorm) 
colnames(a0s.unnorm.post) <-  paste("dataset_", 1:D, sep = "")
a0.unnorm.dt <- melt(a0s.unnorm.post, variable.name = "data_set", value.name = "a0")
a0.unnorm.dt$normalisation <- rep("none", length(a0.unnorm.dt))

a0s.approx.post <- data.frame(a0 = a0.approx) 
colnames(a0s.approx.post) <-  paste("dataset_", 1:D, sep = "")
a0.approx.dt <- melt(a0s.approx.post, variable.name = "data_set", value.name = "a0")
a0.approx.dt$normalisation <- rep("approximate", length(a0.approx.dt))

a0.dt <- rbind(a0.unnorm.dt, a0.approx.dt)

library(ggplot2)

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = data_set, colour = data_set)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, eta[1], nu[1]),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Logistic regression") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist

# ggsave("figure_name.pdf", a0_dist)
###
unnorm.pars <- as.data.frame(
  cbind(extract(unnorm.posterior.lgr, 'alpha')$alpha,
        extract(unnorm.posterior.lgr, 'beta')$beta)
)
names(unnorm.pars) <- c("alpha", paste("beta[", 1:4, "]", sep = ""))

##

approx.norm.pars <- as.data.frame(
  cbind(extract(approx.norm.posterior.lgr, 'alpha')$alpha,
        extract(approx.norm.posterior.lgr, 'beta')$beta)
)
names(approx.norm.pars) <- c("alpha", paste("beta[", 1:4, "]", sep = ""))

##

unnorm.dt <- melt(unnorm.pars, variable.name = "parameter", value.name = "sample")
unnorm.dt$normalisation <- "none"

approx.norm.dt <- melt(approx.norm.pars, variable.name = "parameter", value.name = "sample")
approx.norm.dt$normalisation <- "approximate"

posterior.dt <- rbind(unnorm.dt, approx.norm.dt)

###

true.pars <- data.frame(parameter = names(unnorm.pars), value = c(true.alpha, true.betas))

parameter_posteriors <- ggplot(data = posterior.dt, aes(x = sample, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  scale_x_continuous("") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  facet_wrap(.~parameter,  scales = "free", labeller = label_parsed) +
  geom_vline(data = true.pars, aes(xintercept = value), linetype = "dashed") +
  theme_bw(base_size = 16)

parameter_posteriors

# ggsave("figure_name.pdf", parameter_posteriors)
