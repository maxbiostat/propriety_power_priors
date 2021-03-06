library(npowerPrioR)
source("data_Poisson.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"
prior <- stan_model("stan/simple_Poisson_prior.stan")

po.data <- list(
  N0 = N_0,
  y0 = y_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = NA
)

####################
J <- 20
maxA <- 1
epsilon <- 0.05

adaptive.time <- system.time(
  adaptive.ca0.estimates <- build_grid(compiled.model.prior = prior,
                                       eps = epsilon, M = maxA, J = J, v1 = 10, v2 = 10,
                                       stan.list = po.data, pars = "lambda")
)

## Exporting 
write.csv(adaptive.ca0.estimates$result, row.names = FALSE,
          file = paste("../data/constant_data/Poisson_logCA0_adaptive_J=", J, ".csv", sep = ""))

