source("power_priors_aux.r")

scenario <- 4

source(paste("analyses_Bernoulli/data_Bernoulli_scenario_", scenario, ".r", sep = ""))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_Bernoulli_prior.stan") 

bb.data <- list(
  N0 = N_0,
  y0 = y_0,
  c = cc,
  d = dd,
  a_0 = NA
)

#####################

epsilon <- 0.05
J <- 20
maxA <- 1

source("grid_builder.r")

adaptive.time <- system.time(
  adaptive.ca0.estimates <- build_grid(eps = epsilon,
                                       M = maxA, J = J, v1 = 10, v2 = 10, stan.list = bb.data, pars = "theta")
)
write.csv(adaptive.ca0.estimates$result,
          file = paste("../data/constant_data/Bernoulli_logCA0_scenario_", scenario,
                       "_J=", J, ".csv", sep = ""), row.names = FALSE)


summaries.adaptive <- adaptive.ca0.estimates$summaries
for (j in 1:length(summaries.adaptive)){
  summaries.adaptive[[j]] <- data.frame(summaries.adaptive[[j]], a_0 = adaptive.ca0.estimates$result$a0[j])
} 
summaries.adaptive.dt <- do.call(rbind, summaries.adaptive)


write.csv(summaries.adaptive.dt,
          file = paste("../data/sensitivity_data/priorSensitivity_simpleBernoulli_adaptive_scenario_", scenario,
          "_J=", J, ".csv", sep = ""))