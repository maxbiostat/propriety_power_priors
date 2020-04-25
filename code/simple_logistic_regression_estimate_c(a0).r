scenario <- "B"

source(paste("data_logistic_regression_scenario_", scenario, ".r", sep = ""))
source("power_priors_aux.r")
summary(glm(y_0 ~  X_0, family = "binomial"))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"

compiled.model.prior <- stan_model("stan/simple_logistic_regression_prior.stan")

lgr.data <- list(
  N0 = N_0,
  X0 = X_0,
  P = P,
  y0 = y_0,
  a_0 = NULL
)


###
source("grid_builder.r")
J <- 20
maxA <- 1
epsilon <- 0.05

adaptive.time <- system.time(
  adaptive.ca0.estimates <- build_grid(eps = epsilon,
                                       M = maxA, J = J, v1 = 10, v2 = 10, stan.list = lgr.data, pars = c("alpha", "beta"))
)

write.csv(adaptive.ca0.estimates$result,
          file = paste("../data/constant_data/RegressionLogistic_logCA0_J=", J, "_scenario_", scenario, ".csv", sep = ""), row.names = FALSE)

save(adaptive.ca0.estimates, file = paste("../data/sensitivity_data/RegressionLogistic_logCA0_J=", J, "_scenario_", scenario, ".RData", sep = ""))