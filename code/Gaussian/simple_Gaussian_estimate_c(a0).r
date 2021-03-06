library(npowerPrioR)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

source("data_Gaussian.r")


#### Sampling from the "prior"
prior <- stan_model("stan/simple_Gaussian_prior.stan")

gs.data <- list(
  N0 = N_0,
  y0 = y_0,
  mu0 = mu_0,
  kappa0 = kappa_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = 1
)

####################
J <- 20
maxA <- 10
epsilon <- 0.05

#
uniform.a0s <- seq(epsilon, maxA, length.out = J)
uniform.time <- system.time(
  uniform.ca0.estimates <- create_lc_df(a0_grid = uniform.a0s, compiled.model.prior = prior,
                                        stan.list = gs.data, pars = c("mu", "sigma")) 
)

adaptive.time <- system.time(
  adaptive.ca0.estimates <- build_grid(compiled.model.prior = prior, eps = epsilon,
                                       M = maxA, J = J, v1 = 10, v2 = 10, stan.list = gs.data, pars = c("mu", "sigma"))
)

adaptive.derivOnly.time <- system.time(
  adaptive.ca0.estimates.derivOnly <- build_grid_derivOnly(compiled.model.prior = prior, eps = epsilon,
                                                           M = maxA, J = J, v1 = 10, v2 = 10,
                                                           stan.list = gs.data, pars = c("mu", "sigma"))
)
## Exporting 
write.csv(uniform.ca0.estimates$result, row.names = FALSE,
          file =  paste("../data/constant_data/Gaussian_logCA0_uniform", "_J=", J, ".csv", sep = ""))
write.csv(adaptive.ca0.estimates$result, row.names = FALSE,
          file = paste("../data/constant_data/Gaussian_logCA0_adaptive", "_J=", J, ".csv", sep = ""))
write.csv(adaptive.ca0.estimates.derivOnly$result, row.names = FALSE,
          file = paste("../data/constant_data/Gaussian_logCA0_adaptive_derivOnly", "_J=", J, ".csv", sep = ""))

## Exporting sensitivity stuff

summaries.uniform <- uniform.ca0.estimates$summaries
for (i in 1:length(summaries.uniform)){
  summaries.uniform[[i]] <- data.frame(summaries.uniform[[i]], a_0 = uniform.ca0.estimates$result$a0[i])
} 
summaries.uniform.dt <- do.call(rbind, summaries.uniform)

summaries.adaptive <- adaptive.ca0.estimates$summaries
for (j in 1:length(summaries.adaptive)){
  summaries.adaptive[[j]] <- data.frame(summaries.adaptive[[j]], a_0 = adaptive.ca0.estimates$result$a0[j])
} 
summaries.adaptive.dt <- do.call(rbind, summaries.adaptive)

write.csv(summaries.uniform.dt,
          file = paste("../data/sensitivity_data/priorSensitivity_simpleGaussian_uniform_J=", J, ".csv", sep = ""))

write.csv(summaries.adaptive.dt,
          file = paste("../data/sensitivity_data/priorSensitivity_simpleGaussian_adaptive_J=", J, ".csv", sep = ""))