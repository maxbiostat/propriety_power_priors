source("data_regression_NIG.r")
summary(lm(y_0 ~ -1 + X_0))

library(rstan)
rstan_options(auto_write = TRUE)


lm.data <- list(
  N0 = N_0,
  X0 = X_0,
  y0 = y_0,
  P = ncol(X_0),
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  a_0 = .15
)

invlambda0 <- solve(lm.data$lambda_0)


############
### Unnormalised posterior

compiled.model.unnorm.posterior <- stan_model("stan/simple_linear_regression_NIG_posterior_unnormalised.stan")

lm.data.forposterior <- list(
  N0 = N_0,
  X0 = X_0,
  y0 = y_0,
  P = ncol(X_0),
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  eta = eta,
  nu = nu,
  X = X,
  N = N,
  y = y
)

unnorm.posterior.lm <- sampling(compiled.model.unnorm.posterior, data = lm.data.forposterior)
unnorm.posterior.lm
pairs(unnorm.posterior.lm, pars = c("a_0", "beta"))

### The approximately normalised prior

J <- 20

constant_data <- read.csv(paste("../data/constant_data/RegressionNIG_logCA0_adaptive_J=", J, ".csv", sep = ""))

library(mgcv)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

# plot(constant_data$a0, predict.gam(fit.gam), type = "l")
# points(lc_a0 ~ a0 , constant_data)

approx.normalised.model <- stan_model("stan/simple_linear_regression_NIG_posterior_normalised_approximate.stan")

Ks <- c(50, 100, 1000, 10000)

get_posterior_K <- function(K){
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  lm.data.forposterior$pred_grid_x <- a0_grid$a0
  lm.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.lm <- sampling(approx.normalised.model,
                                          data = lm.data.forposterior, refresh = 0)
  return(approx.norm.posterior.lm)
}

system.time(
  all.approximates <- lapply(Ks, get_posterior_K)  
)

## Now the normalised version

normalised.model <- stan_model("stan/simple_linear_regression_NIG_posterior_normalised.stan")

norm.posterior.lm <- sampling(normalised.model, data = lm.data.forposterior)
norm.posterior.lm
pairs(norm.posterior.lm, pars = c("a_0", "beta"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.lm, 'a_0')$a_0

a0.approx.list <- lapply(Ks, function(k){
  i <- match(k, Ks)
  a0 <- extract(all.approximates[[i]], 'a_0')$a_0
  return(data.frame(a0 = a0, normalisation = paste("K=", k, sep = "")))
} )

a0.norm <- extract(norm.posterior.lm, 'a_0')$a_0

a0.dt <- 
  rbind(
    data.frame(a0 = c(a0.unnorm, a0.norm),
               normalisation = c( rep("none", length(a0.unnorm)),
                                  rep("exact", length(a0.norm)) )
    ),
    do.call(rbind, a0.approx.list)
  )

library(ggplot2)

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density(alpha = .4) +
  stat_function(fun = function(x) dbeta(x, eta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Linear regression") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist

ggsave(filename = "../figures/a0_posterior_RegressionNIG.pdf", plot = a0_dist, dpi = 300)

a0_dist_norm <- ggplot(subset(a0.dt, normalisation != "none"),
                       aes(x = normalisation, y = a0,
                           fill = normalisation, colour = normalisation)) +
  geom_boxplot(alpha = .4) +
  ggtitle("Linear regression") +
  scale_y_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist_norm

ggsave("../figures/a0_posterior_RegressionNIG_normalisation_comparison.pdf", a0_dist_norm)

###
unnorm.pars <- as.data.frame(
  cbind(extract(unnorm.posterior.lm, 'beta')$beta,
        extract(unnorm.posterior.lm, 'sigma_sq')$sigma_sq)
)
names(unnorm.pars) <- c(paste("beta[", 1:4, "]", sep = ""), "sigma^2")

##

app.norm.pars.list <- lapply(all.approximates, function(x) {
  pars <- as.data.frame(
    cbind(extract(x, 'beta')$beta,
          extract(x, 'sigma_sq')$sigma_sq)
  )
  names(pars) <- c(paste("beta[", 1:4, "]", sep = ""), "sigma^2")
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
norm.pars <- as.data.frame(
  cbind(extract(norm.posterior.lm, 'beta')$beta,
        extract(norm.posterior.lm, 'sigma_sq')$sigma_sq)
)
names(norm.pars) <- c(paste("beta[", 1:4, "]", sep = ""), "sigma^2")

unnorm.dt <- melt(unnorm.pars, variable.name = "parameter", value.name = "sample")
unnorm.dt$normalisation <- "none"
norm.dt <- melt(norm.pars, variable.name = "parameter", value.name = "sample")
norm.dt$normalisation <- "exact"

posterior.dt <- rbind(unnorm.dt, norm.dt, do.call(rbind, app.norm.pars.df.list))
posterior.dt$normalisation <- factor(posterior.dt$normalisation,
                                     levels = c("none", paste("K=", Ks, sep = ""), "exact") )

###

true.pars <- data.frame(parameter = names(norm.pars), value = c(true.betas, sy^2))
  
parameter_posteriors <- ggplot(data = posterior.dt, aes(x = sample, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  scale_x_continuous("") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  facet_wrap(.~parameter,  scales = "free", labeller = label_parsed) +
  geom_vline(data = true.pars, aes(xintercept = value), linetype = "dashed") +
  theme_bw(base_size = 16)

parameter_posteriors

ggsave("../figures/parameter_posterior_Regression.pdf", parameter_posteriors)
