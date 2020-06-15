source("data_logistic_regression_multiple_historical.r")
source("../power_priors_aux.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"

compiled.model.prior <- stan_model("stan/simple_logistic_regression_prior.stan")
data.lists <- vector(D, mode = "list")
exclude_last <- function(dt) dt[, -ncol(dt)]
get_last <- function(dt) dt[, ncol(dt)]
  
for (d in 1:D){
  data.lists[[d]] <- list(
    N0 = N_0,
    X0 = exclude_last(dataSets[[d]]),
    P = P,
    y0 = get_last(dataSets[[d]]),
    a_0 = NULL
  )
}

save(data.lists, "data/data_lists.RData")

###
source("../grid_builder.r")
J <- 20
maxA <- 1
epsilon <- 0.05
Alls <- vector(D, mode = "list")

build.time <- system.time( ## could be made faster by parallelising at the dataset level (using parallel::mclapply() for instance)
  for(d in 1:D){
    Alls[[d]] <- build_grid(eps = epsilon, M = maxA, J = J, v1 = 10, v2 = 10,
                            stan.list = data.lists[[d]], pars = c("alpha", "beta"))
  }
)
build.time

save(Alls, file = paste("data/RegressionLogistic_multiple_historical_logCA0_J=", J, ".RData", sep = ""))

for (d in 1:D){
  res <- Alls[[d]]$result
  res$dataset <- paste("dataset_", d, sep = "")
  write.csv(res,
            file = paste("data/RegressionLogistic_logCA0_J=", J, "_dataset_", d, ".csv", sep = ""), row.names = FALSE)
}

