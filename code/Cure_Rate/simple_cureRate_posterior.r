source("data_cure_rate_real.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

############
### Unnormalised posterior

compiled.model.unnorm.posterior <- stan_model("stan/cure_rate_posterior_unnormalised.stan")

cr.data.forposterior <- list(
  N0 = N_0,
  P = ncol(X_0),
  X0 = X_0,
  y0 = Y_0_cens,
  Delta0 = Delta_0,
  N = N,
  X = X,
  y = Y_cens,
  Delta = Delta,
  mu_0 = mu0,
  sigma_0 = s0,
  delta_0 = d0,
  tau_0 = tau0,
  sigma_beta = 10,
  eta = eta,
  nu = nu
)

unnorm.posterior.cr <- sampling(compiled.model.unnorm.posterior, data = cr.data.forposterior)
unnorm.posterior.cr
pairs(unnorm.posterior.cr, pars = c("a_0", "beta"))
pairs(unnorm.posterior.cr, pars = c("a_0", "alpha", "lambda"))

### The approximately normalised prior
 J <- 20

constant_data <- read.csv("../data/constant_data/CureRateRealData_logCA0_adaptive.csv")
library(mgcv)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)


approx.normalised.model <- stan_model("stan/cure_rate_posterior_normalised_approximate.stan")

Ks <- c(100, 1000, 10000, 20000)

get_posterior_K <- function(K){
  cat("Doing k=", K, "\n")
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  cr.data.forposterior$pred_grid_x <- a0_grid$a0
  cr.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.cr <- sampling(approx.normalised.model,
                                          data = cr.data.forposterior,
                                       control = list(adapt_delta = .95, max_treedepth = 15),
                                       refresh = 0)
  return(approx.norm.posterior.cr)
}

system.time(
  all.approximates <- lapply(Ks, get_posterior_K)  
)

all.approximates[[4]]

pairs(all.approximates[[4]], pars = c("a_0", "beta"))
pairs(all.approximates[[4]], pars = c("a_0", "alpha", "lambda"))


##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.cr, 'a_0')$a_0

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
  stat_function(fun = function(x) dbeta(x, eta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Logistic regression") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist

ggsave("../figures/a0_posterior_CureRate_prior2.pdf", a0_dist)

# subset(a0.dt, normalisation != "unnormalised")
a0_dist_norm <- ggplot(a0.dt,
                       aes(x = normalisation, y = a0,
                           fill = normalisation, colour = normalisation)) +
  geom_boxplot(alpha = .4) +
  ggtitle("Logistic regression") +
  scale_y_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist_norm

ggsave("../figures/a0_posterior_CureRate_normalisation_comparison_prior2.pdf", a0_dist_norm)

###
unnorm.pars <- as.data.frame(
  cbind(extract(unnorm.posterior.cr, 'alpha')$alpha,
        extract(unnorm.posterior.cr, 'lambda')$lambda,
  extract(unnorm.posterior.cr, 'beta')$beta)
)
names(unnorm.pars) <- c("alpha", "lambda", paste("beta[", 1:4, "]", sep = ""))

##

app.norm.pars.list <- lapply(all.approximates, function(x) {
  pars <- as.data.frame(
    cbind(extract(x, 'alpha')$alpha,
          extract(x, 'lambda')$lambda,
          extract(x, 'beta')$beta
          )
  )
  names(pars) <- c("alpha", "lambda", paste("beta[", 1:4, "]", sep = ""))
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
unnorm.dt$normalisation <- "unnormalised"

posterior.dt <- rbind(unnorm.dt, do.call(rbind, app.norm.pars.df.list))
posterior.dt$normalisation <- factor(posterior.dt$normalisation,
                                     levels = c("unnormalised", paste("K=", Ks, sep = "")) )

###
parameter_posteriors <- ggplot(data = posterior.dt, aes(x = sample, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  scale_x_continuous("") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  facet_wrap(.~parameter,  scales = "free", labeller = label_parsed) +
  theme_bw(base_size = 16)

parameter_posteriors

ggsave("../figures/parameter_posterior_CureRate_prior2.pdf", parameter_posteriors)
