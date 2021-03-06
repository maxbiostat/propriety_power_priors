source("power_priors_aux.r")
source("data_cure_rate_real.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/cure_rate_prior.stan")

cr.data <- list(
  n = N_0,
  p = ncol(X_0),
  y = Y_0_cens,
  X = X_0,
  delta = Delta_0,
  mu_0 = mu0,
  sigma_0 = s0,
  delta_0 = d0,
  tau_0 = tau0,
  sigma_beta = 10,
  a_0 = NA
)
####################
J <- 20
maxA <- 1
epsilon <- 0.05

source("grid_builder.r")
#
adaptive.time <- system.time(
  adaptive.ca0.estimates <- build_grid(eps = epsilon,
                                       M = maxA, J = J, v1 = 10, v2 = 10, stan.list = cr.data,
                                       pars = c("alpha", "lambda", "beta"), strict = TRUE)
)
write.csv(adaptive.ca0.estimates$result, row.names = FALSE,
          file = "../data/constant_data/CureRateRealData_logCA0_adaptive.csv")
save(adaptive.ca0.estimates, file = "../data/sensitivity_data/cure_rate_real.RData")