### In this script we simply fit the cure rate model to the E1684 study data to look at parameter estimates  
### We'll fit the model and compare estimates with Chen, Ibrahim and Sinha (1999) [CIS99]
historical.data <- subset(read.csv("../data/e1684_and_e1690_data.csv", sep = "\t"),
                          study  == "1684" & survtime > 0)
summary(historical.data)

stand_age <- (historical.data$age -mean(historical.data$age))/sd(historical.data$age)
N_0 <- length(stand_age)
X_0 <- cbind(rep(1, N_0), # x0 - dummy for intercept
             stand_age, # x1 - (standardised) age 
             historical.data$sex, # x2 - gender (sex)
             historical.data$perform # x3 performance status
             )
summary(X_0)
Y_0_cens <- historical.data$survtime
Delta_0 <- historical.data$scens
## Hyperparameters
mu0 <- 0
s0 <- 10

d0 <- 1
tau0 <- 1

########################3
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling 
simple.model.cov <- stan_model("stan/cure_rate.stan")

E1684.data <- list(
  n = N_0,
  p = ncol(X_0),
  y = Y_0_cens,
  X = X_0,
  delta = Delta_0,
  mu_0 = mu0,
  sigma_0 = s0,
  delta_0 = d0,
  tau_0 = tau0,
  sigma_beta = 10
)

optimizing(simple.model.cov, data = E1684.data, verbose = TRUE)

mcmc <- sampling(simple.model.cov, data = E1684.data,
                 control = list(adapt_delta = 0.99, max_treedepth = 15, metric = "dense_e"))
mcmc

pairs(mcmc, pars = c("lambda", "alpha", "beta"))

# pars <- extract(mcmc, pars = c("lambda", "alpha", "beta"))
# M <- nrow(pars$beta)
# 
# tcens <- 11
# source("cure_rate_aux.r")
# y_rep <- matrix(NA, ncol = E1684.data$n, nrow = M)
# for(i in 1:M){
#   y_rep[i, ] <- generate_n_data_points(n = E1684.data$n,
#                                        beta = pars$beta[i, ],
#                                        alpha = pars$alpha[i],
#                                        lambda= pars$lambda[i],
#                                        X = E1684.data$X,
#                                        tcs = tcens)$y_censored
# }
# 
# library(bayesplot)
# 
# plot0 <- ppc_ecdf_overlay(E1684.data$y, y_rep[1:min(100, round(M * .2)), ], discrete = TRUE) +
#   # ggtitle("")  +
#   NULL
# plot0
# 
# plot1 <- ppc_stat(E1684.data$y, y_rep[1:M, ], stat = "mean", binwidth = 0.005) +
#   # ggtitle("")  +
#   NULL
# plot1
# 
# plot2 <- ppc_stat(E1684.data$y, y_rep[1:M, ], stat = "sd", binwidth = 0.005) +
#   # ggtitle("")  +
#   NULL
# plot2
# 
# q90 <- function(x) quantile(x, probs = .90)
# plot3 <- ppc_stat(E1684.data$y, y_rep[1:M, ], stat = "q90", binwidth = 0.005) +
#   # ggtitle("")  +
#   NULL
# plot3
