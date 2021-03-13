######## Simple Gaussian distribution example.
######## This script will illustrate the effect of (not) including the normalising constant c(a0).
######## It will also show the effect of various approximations to c(a0)

source("data_Gaussian.r")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#####################
####### Posterior using unnormalised power prior
#####################

## New data
### will draw from the same process for simplicity

unnormalised.model <- stan_model("stan/simple_Gaussian_posterior_unnormalised.stan")

gs.data.forposterior <- list(
  N0 = N_0,
  y0 = y_0,
  mu0 = mu_0,
  kappa0 = kappa_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  eta = eta,
  nu = nu,
  N = N,
  y = y
)

get_pars_gaussian <- function(y0, n0, alpha0, beta0, m0, k0, a_0){
  nstar <- a_0 * n0
  ybar <- mean(y0)
  s <- mean( (y0-ybar)^2 )
  kappa_n <- k0 + nstar  
  alpha_n <- alpha0 + nstar/2
  beta_n <- beta0 + .5 * (nstar * s +  (k0 * nstar * (ybar - m0)^2 )/kappa_n ) 
  return(list(
    alpha_n = alpha_n,
    beta_n = beta_n,
    mu_n = (k0 * m0  + a_0 * n0 * ybar)/(k0 + a_0 * n0)
  ))
}

post_pars <- get_pars_gaussian(y0 = y_0, n0 = N_0, alpha0 = alpha_0, beta0 = beta_0, m0 = mu_0, k0 = kappa_0,
                  a_0 = 1)

unnorm.posterior.gaussian <- sampling(unnormalised.model, data = gs.data.forposterior)
unnorm.posterior.gaussian
pairs(unnorm.posterior.gaussian, pars = c("a_0", "mu", "tau"))

### The approximately normalised prior

J <- 20

constant_data <- read.csv(paste("../../data/constant_data/Gaussian_logCA0_adaptive_J=", J, ".csv", sep = ""))
library(mgcv)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)


approx.normalised.model <- stan_model("stan/simple_Gaussian_posterior_normalised_approximate.stan")

Ks <- c(50, 100, 1000, 10000)

get_posterior_K <- function(K){
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))
  gs.data.forposterior$pred_grid_x <- a0_grid$a0
  gs.data.forposterior$pred_grid_y <- a0_grid$lc_pred
  approx.norm.posterior.gauss <- sampling(approx.normalised.model,
                                          data = gs.data.forposterior, refresh = 0)
  return(approx.norm.posterior.gauss)
}

system.time(
  all.approximates <- lapply(Ks, get_posterior_K)  
)

### Now the normalised prior
normalised.model <- stan_model("stan/simple_Gaussian_posterior_normalised.stan")

norm.posterior.gauss <- sampling(normalised.model, data = gs.data.forposterior)
norm.posterior.gauss
pairs(norm.posterior.gauss, pars = c("a_0", "mu", "tau"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.gaussian, 'a_0')$a_0

a0.approx.list <- lapply(Ks, function(k){
  i <- match(k, Ks)
  a0 <- extract(all.approximates[[i]], 'a_0')$a_0
  return(data.frame(a0 = a0, normalisation = paste("K=", k, sep = "")))
} )

a0.norm <- extract(norm.posterior.gauss, 'a_0')$a_0

a0.dt <- 
  rbind(
    data.frame(a0 = c(a0.unnorm, a0.norm),
               normalisation = c( rep("none", length(a0.unnorm)),
                                  rep("exact", length(a0.norm)) )
    ),
    do.call(rbind, a0.approx.list)
  )

a0.dt$normalisation <- factor(a0.dt$normalisation,
                              levels = c("none",  paste("K=", Ks, sep = ""), "exact") )

library(ggplot2)

if(copy){
  get_l_a0_gaussian <- function(y0, n0, alpha0, beta0, m0, k0, a_0){
    nstar <- a_0 * n0
    ybar <- mean(y0)
    s <- mean( (y0-ybar)^2 )
    kappa_n <- k0 + nstar  
    alpha_n <- alpha0 + nstar/2
    beta_n <- beta0 + .5 * (nstar * s +  (k0 * nstar * (ybar - m0)^2 )/kappa_n ) 
    ans <- lgamma(alpha_n)-lgamma(alpha0)
    ans <- ans + alpha0 * log(beta0) - alpha_n * log(beta_n)
    ans <- ans + .5 *( log(k0) - log(kappa_n) )-nstar/2 * log(2*pi)
    return(ans)
  }
  ############
  l_a0 <- function(x) {
    get_l_a0_gaussian(
      y0 = gs.data.forposterior$y0,
      n0 = gs.data.forposterior$N0,
      alpha0 = gs.data.forposterior$alpha0,
      beta0 = gs.data.forposterior$beta0,
      m0 = gs.data.forposterior$mu0,
      k0 = gs.data.forposterior$kappa0,
      a_0 = x
    )
  }
  l_a0 <- Vectorize(l_a0)
  
  special_marginal_a0 <- function(x, eta = 1, nu = 1, log = FALSE){
    ans <- dbeta(x, shape1 = eta, shape2 = nu, log = TRUE) + l_a0(x + 1) - l_a0(x)
    if(!log) ans <- exp(ans)
    return(ans)
  }
  K <- integrate(special_marginal_a0, 0, 1)$value
  special_marginal_a0_norm <- function(x) special_marginal_a0(x)/K
  special_marginal_a0_norm <- Vectorize(special_marginal_a0_norm)
}

a0_dist <- ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, eta, nu),
                geom = "line", colour = "black", linetype = "longdash") + 
  ggtitle("Gaussian") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom",
        legend.justification = "centre",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0))


if(copy) a0_dist <- a0_dist +   stat_function(fun = special_marginal_a0_norm,
                                              geom = "line", colour = "black", linetype = "solid") 

a0_dist

ggsave(paste("../../figures/a0_posterior_Gaussian_J=", J, ".pdf", sep = ""), a0_dist, dpi = 300)

a0_dist_norm <- ggplot(subset(a0.dt, normalisation != "none"),
                       aes(x = normalisation, y = a0,
                           fill = normalisation, colour = normalisation)) +
  geom_boxplot(alpha = .4) +
  ggtitle("Gaussian") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom",
        legend.justification = "centre",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0))


a0_dist_norm

ggsave(paste("../../figures/a0_posterior_Gaussian_normalisation_comparison_J=", J, ".pdf", sep = ""),
       a0_dist_norm, dpi = 300)

###
unnorm.pars <- as.data.frame(
  cbind(extract(unnorm.posterior.gaussian, 'mu')$mu,
        extract(unnorm.posterior.gaussian, 'tau')$tau)
)
names(unnorm.pars) <- c("mu", "tau")

app.norm.pars.list <- lapply(all.approximates, function(x) {
  pars <- as.data.frame(
    cbind(extract(x, 'mu')$mu,
          extract(x, 'tau')$tau)
  )
  names(pars) <- c("mu", "tau")
  return(pars)
})

library(reshape2)
app.norm.pars.df.list <- lapply(seq_along(app.norm.pars.list), function(i) {
  y <- app.norm.pars.list[[i]]
  dt <- melt(y, variable.name = "parameter", value.name = "sample")
  dt$normalisation <- paste("K=", Ks[i] ,sep = "")
  return(dt)
} )

#
norm.pars <- as.data.frame(
  cbind(extract(norm.posterior.gauss, 'mu')$mu,
        extract(norm.posterior.gauss, 'tau')$tau)
)
names(norm.pars) <- c("mu", "tau")

unnorm.dt <- melt(unnorm.pars, variable.name = "parameter", value.name = "sample")
unnorm.dt$normalisation <- "none"
norm.dt <- melt(norm.pars, variable.name = "parameter", value.name = "sample")
norm.dt$normalisation <- "exact"
posterior.dt <- rbind(unnorm.dt, norm.dt, do.call(rbind, app.norm.pars.df.list))

posterior.dt$normalisation <- factor(posterior.dt$normalisation,
                                     levels = c("none",  paste("K=", Ks, sep = ""), "exact") )

###

post.pars <- get_pars_gaussian(y0 = y_0, n0 = N_0, alpha0 = alpha_0,
                               beta0 = beta_0, m0 = mu_0, k0 = kappa_0,
                               a_0 = .95)

post.means <- data.frame(parameter = names(norm.pars), value = c(post.pars$mu_n, post.pars$alpha_n/post.pars$beta_n))

parameter_posteriors <- ggplot(data = posterior.dt,
                               aes(x = sample, colour = normalisation, fill = normalisation)) +
  geom_density(alpha = .4) +
  scale_x_continuous("") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  facet_wrap(.~parameter,  scales = "free", labeller = label_parsed) +
  # geom_vline(data = post.means, aes(xintercept = value), linetype = "dashed") +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom",
        legend.justification = "centre",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0))


parameter_posteriors

ggsave(filename = paste("../../figures/parameter_posterior_Gaussian_J=", J, ".pdf", sep = ""),
       plot = parameter_posteriors, dpi = 300)
