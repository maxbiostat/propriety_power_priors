source("power_priors_aux.r")
source("data_Gaussian.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_Gaussian_posterior_unnormalised_fixed_a0.stan") ## lying a bit, this is the posterior

gs.data <- list(
  N0 = N_0,
  y0 = y_0,
  mu0 = mu_0,
  kappa0 = kappa_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  N = N,
  y = y
)
####################
J <- 20
maxA <- 10
epsilon <- 0.05

source("grid_builder.r")
#
uniform.a0s <- seq(epsilon, maxA, length.out = J)
uniform.time <- system.time(
  uniform.ca0.estimates <- create_lc_df(a0_grid = uniform.a0s, stan.list = gs.data, pars = c("mu", "sigma")) 
)
## for the sensitivity analysis all we need are estimates of parameters, no need for l(a_0)
adaptive.derivOnly.time <- system.time(
    adaptive.ca0.estimates.derivOnly <- build_grid_derivOnly(eps = epsilon,
                                                           M = maxA, J = J, v1 = 10, v2 = 10,
                                                           stan.list = gs.data, pars = c("mu", "sigma"))
)
## Exporting 
write.csv(uniform.ca0.estimates$result, row.names = FALSE,
          file =  paste("../data/constant_data/Gaussian_posterior_logCA0_uniform", "_J=", J,
                        ".csv", sep = ""))
write.csv(adaptive.ca0.estimates.derivOnly$result, row.names = FALSE,
          file = paste("../data/constant_data/Gaussian_posterior_logCA0_adaptive_derivOnly", "_J=", J,
                       ".csv", sep = ""))

## Exporting sensitivity stuff

summaries.uniform <- uniform.ca0.estimates$summaries
for (i in 1:length(summaries.uniform)){
  summaries.uniform[[i]] <- data.frame(summaries.uniform[[i]], a_0 = uniform.ca0.estimates$result$a0[i])
} 
summaries.uniform.dt <- do.call(rbind, summaries.uniform)

summaries.adaptive <- adaptive.ca0.estimates.derivOnly$summaries
for (j in 1:length(summaries.adaptive)){
  summaries.adaptive[[j]] <- data.frame(summaries.adaptive[[j]],
                                        a_0 = adaptive.ca0.estimates.derivOnly$result$a0[j])
} 
summaries.adaptive.dt <- do.call(rbind, summaries.adaptive)

write.csv(summaries.uniform.dt,
          file = paste("../data/sensitivity_data/posteriorSensitivity_simpleGaussian_uniform_J=", J, ".csv", sep = ""))

write.csv(summaries.adaptive.dt,
          file = paste("../data/sensitivity_data/posteriorSensitivity_simpleGaussian_adaptive_J=", J, ".csv", sep = ""))