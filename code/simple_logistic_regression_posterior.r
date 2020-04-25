scenario <- "B"

source(paste("data_logistic_regression_scenario_", scenario, ".r", sep = ""))
source("power_priors_aux.r")

summary(glm(y_0 ~ X_0, family = "binomial"))
summary(glm(y ~ X, family = "binomial"))
summary(glm(c(y_0, y) ~ rbind(X_0, X), family = "binomial"))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

############
### Unnormalised posterior

compiled.model.unnorm.posterior <- stan_model("stan/simple_logistic_regression_posterior_unnormalised.stan")

lgr.data.forposterior <- list(
  N0 = N_0,
  X0 = X_0,
  P = P,
  y0 = y_0,
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

J <- 20
constant_data <- read.csv(paste("../data/constant_data/RegressionLogistic_logCA0_J=", J, "_scenario_", scenario, ".csv", sep = ""))

library(mgcv)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)


approx.normalised.model <- stan_model("stan/simple_logistic_regression_posterior_approximate.stan")

Ks <- c(100, 1000, 10000, 20000)

get_posterior_K <- function(K){
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  lgr.data.forposterior$pred_grid_x <- a0_grid$a0
  lgr.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.lgr <- sampling(approx.normalised.model,
                                          data = lgr.data.forposterior, refresh = 0)
  return(approx.norm.posterior.lgr)
}

system.time(
  all.approximates <- lapply(Ks, get_posterior_K)  
)


##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.lgr, 'a_0')$a_0

a0.approx.list <- lapply(Ks, function(k){
  i <- match(k, Ks)
  a0 <- extract(all.approximates[[i]], 'a_0')$a_0
  return(data.frame(a0 = a0, normalisation = paste("K=", k, sep = "")))
} )

a0.dt <- 
  rbind(
    data.frame(a0 = a0.unnorm,
               normalisation = rep("none", length(a0.unnorm)) 
    ),
    do.call(rbind, a0.approx.list)
  )

library(ggplot2)

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, eta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Logistic regression") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist

ggsave(paste("../figures/a0_posterior_RegressionLogistic_J=", J, "_scenario_", scenario, ".pdf", sep = ""), a0_dist)

a0_dist_norm <- ggplot(subset(a0.dt, normalisation != "none"),
                       aes(x = normalisation, y = a0,
                           fill = normalisation, colour = normalisation)) +
  geom_boxplot(alpha = .4) +
  ggtitle("Logistic regression") +
  scale_y_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist_norm

ggsave(paste("../figures/a0_posterior_RegressionLogistic_normalisation_comparison_J=", J, "_scenario_", scenario, ".pdf", sep = ""), a0_dist_norm)

###
unnorm.pars <- as.data.frame(
  cbind(extract(unnorm.posterior.lgr, 'alpha')$alpha,
  extract(unnorm.posterior.lgr, 'beta')$beta)
)
names(unnorm.pars) <- c("alpha", paste("beta[", 1:4, "]", sep = ""))

##

app.norm.pars.list <- lapply(all.approximates, function(x) {
  pars <- as.data.frame(
    cbind(extract(x, 'alpha')$alpha,
          extract(x, 'beta')$beta
          )
  )
  names(pars) <- c("alpha", paste("beta[", 1:4, "]", sep = ""))
  return(pars)
})

library(reshape2)
app.norm.pars.df.list <- lapply(seq_along(app.norm.pars.list), function(i) {
  y <- app.norm.pars.list[[i]]
  dt <- melt(y, variable.name = "parameter", value.name = "sample")
  dt$normalisation <- paste("K=", Ks[i] ,sep = "")
  return(dt)
} )

##

unnorm.dt <- melt(unnorm.pars, variable.name = "parameter", value.name = "sample")
unnorm.dt$normalisation <- "none"

posterior.dt <- rbind(unnorm.dt, do.call(rbind, app.norm.pars.df.list))
posterior.dt$normalisation <- factor(posterior.dt$normalisation,
                                     levels = c("none", paste("K=", Ks, sep = "")) )

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

ggsave(paste("../figures/parameter_posterior_LogisticRegression_J=", J, "_scenario_", scenario, ".pdf", sep = ""), parameter_posteriors)
