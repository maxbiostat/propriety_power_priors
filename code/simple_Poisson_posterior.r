source("data_Poisson.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)


constant_data <- read.csv("../data/constant_data/Poisson_logCA0_adaptive_J=20.csv")
maxA <- max(constant_data$a0)

#####################
####### Posterior using unnormalised power prior
#####################

compiled.model.unnorm.posterior <- stan_model("stan/simple_Poisson_unnormalised_posterior.stan")

po.data.forposterior <- list(
  N0 = N_0,
  y0 = y_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  eta = eta,
  nu = nu,
  N = N,
  y = y
)

unnorm.posterior.poisson <- sampling(compiled.model.unnorm.posterior, data = po.data.forposterior)
unnorm.posterior.poisson
pairs(unnorm.posterior.poisson, pars = c("a_0", "lambda"))

### The approximately normalised posterior

library(mgcv)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0), data = constant_data)

approx.normalised.model <- stan_model("stan/simple_Poisson_posterior_normalised_approximate.stan")

Ks <- c(50, 100, 1000, 10000)

get_posterior_K <- function(K){
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  po.data.forposterior$pred_grid_x <- a0_grid$a0
  po.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.poisson <- sampling(approx.normalised.model,
                                            data = po.data.forposterior, refresh = 0)
  return(approx.norm.posterior.poisson)
}

system.time(
  all.approximates <- lapply(Ks, get_posterior_K)  
)

### Now the normalised posterior
compiled.model.norm.posterior <- stan_model("stan/simple_Poisson_normalised_posterior.stan")

norm.posterior.poisson <- sampling(compiled.model.norm.posterior, data = po.data.forposterior)
norm.posterior.poisson
pairs(norm.posterior.poisson, pars = c("a_0", "lambda"))


##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.poisson, 'a_0')$a_0

a0.approx.list <- lapply(Ks, function(k){
  i <- match(k, Ks)
  a0 <- extract(all.approximates[[i]], 'a_0')$a_0
  return(data.frame(a0 = a0, normalisation = paste("K=", k, sep = "")))
} )

a0.norm <- extract(norm.posterior.poisson, 'a_0')$a_0

a0.dt <- 
  rbind(
    data.frame(a0 = c(a0.unnorm, a0.norm),
               normalisation = c( rep("none", length(a0.unnorm)),
                                  rep("exact", length(a0.norm)) )
    ),
    do.call(rbind, a0.approx.list)
  )

a0.dt$normalisation <- factor(a0.dt$normalisation,
                                     levels = c("none", paste("K=", Ks, sep = ""), "exact") )


library(ggplot2)

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, eta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  # ggtitle("Gaussian") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20)

a0_dist


ggsave("../figures/a0_posterior_Poisson.pdf", a0_dist)
###

unnorm.lambda.dt <- data.frame(lambda = extract(unnorm.posterior.poisson, 'lambda')$lambda)
unnorm.lambda.dt$normalisation <- "none"

app.norm.pars.list <- lapply(seq_along(all.approximates), function(i) {
  x <- all.approximates[[i]]
  pars <- data.frame(
    lambda = extract(x, 'lambda')$lambda,
    normalisation = paste("K=", Ks[i], sep = "")
  )
  return(pars)
})

norm.lambda.dt <- data.frame(lambda = extract(norm.posterior.poisson, 'lambda')$lambda)
norm.lambda.dt$normalisation <- "exact"


posterior.dt <- rbind(unnorm.lambda.dt, norm.lambda.dt, do.call(rbind, app.norm.pars.list))
posterior.dt$normalisation <- factor(posterior.dt$normalisation,
                                     levels = c("none", paste("K=", Ks, sep = ""), "exact") )

lambda_posterior <- ggplot(data = posterior.dt, aes(x = lambda, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  # stat_function(fun = function(x) dgamma(x, shape = a_0*sum(y_0) + alpha_0, rate = a_0*N_0 + beta_0),
  #               geom = "line", colour = "black", linetype = "solid") + 
  scale_x_continuous(expression(lambda), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +
  geom_vline(xintercept = true.lambda, linetype = "dashed") +
  theme_bw(base_size = 20)

lambda_posterior

ggsave("../figures/parameter_posterior_Poisson.pdf", lambda_posterior)
