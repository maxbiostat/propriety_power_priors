scenario <- 1

source(paste("data_regression_NIG_scenario_", scenario, ".r", sep = ""))
source("power_priors_aux.r")
summary(lm(y_0 ~ -1 + X_0))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"

compiled.model.prior <- stan_model("stan/simple_linear_regression_NIG_prior.stan")

lm.data <- list(
  N0 = N_0,
  P = P,
  X0 = X_0,
  y0 = y_0,
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  a_0 = NULL
)
#########
J <- 20
maxA <- 1
epsilon <- 0.05

source("grid_builder.r")

adaptive.time <- system.time(
  adaptive.ca0.estimates <- build_grid(eps = epsilon,
                                       M = maxA, J = J, v1 = 10, v2 = 10, stan.list = lm.data, pars = c("beta", "sigma_sq"))
)
## Exporting 
write.csv(adaptive.ca0.estimates$result, row.names = FALSE,
          file = paste("../data/constant_data/RegressionNIG_logCA0_adaptive_J=", J,"_scenario_", scenario, ".csv",  sep = ""))

save(adaptive.ca0.estimates,
     file = paste("../data/constant_data/RegressionNIG_priorParameterEstimates_adaptive_J=", J,"_scenario_", scenario, ".RData",  sep = ""))
