source("power_priors_aux.r")
source("Gaussian_data.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_gaussian_prior.stan")

gs.data <- list(
  N0 = N_0,
  y0 = y_0,
  mu0 = mu_0,
  kappa0 = kappa_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = NA
)
get_c_a0_brute_force <- function(x){
  gs.data <- list(
    N0 = N_0,
    y0 = y_0,
    mu0 = mu_0,
    kappa0 = kappa_0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    a_0 = x
  )
  prior.gauss <- sampling(compiled.model.prior, data = gs.data, refresh = 0)
  res <- list(
    chain = prior.gauss,
    logml = bridgesampling::bridge_sampler(prior.gauss, silent = TRUE)$logml
  )
  return(res)
}
###
get_c_a0_gaussian <- function(y0, n0, alpha0, beta0, m0, k0, a_0){
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
c_a0 <- function(x) {
  get_c_a0_gaussian(
    y0 = gs.data$y0,
    n0 = gs.data$N0,
    alpha0 = gs.data$alpha0,
    beta0 = gs.data$beta0,
    m0 = gs.data$mu0,
    k0 = gs.data$kappa0,
    a_0 = x
  )
}
c_a0 <- Vectorize(c_a0)
####################
J <- 15
maxA <- 10
# est_a0s <- pracma::logseq(0.01, maxA, n = J)
est_a0s <- c(seq(0.05, 1, length.out = J-3), seq(1.2, maxA, length.out = 3))
# est_a0s <- seq(0, maxA, length.out = J)
system.time(
  results <- lapply(est_a0s, get_c_a0_brute_force)  
)
mls <- unlist(lapply(results, function(x) x$logml))
## including the (0, 0) point
a0s <- c(0, est_a0s)
mls <- c(0, mls)

constant_data <- data.frame(a0 = a0s, lc_a0 = mls)
write.csv(constant_data, row.names = FALSE, "../data/Gaussian_logCA0.csv")

#######
K <- 50
pred_a0s <- seq(0, maxA, length.out = K)
degree <- 2

fit.poly <- lm(lc_a0 ~ poly(a0, degree), data = constant_data)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0), data = constant_data)

##########################
preds.poly <- get_band(fit.poly, xpred = pred_a0s)
preds.gam  <- get_band(fit.gam, xpred = pred_a0s)

preds.list <- list(
  polynomial = data.frame (a0 = pred_a0s , lca0 = preds.poly$mean, lwr = preds.poly$lwr, upr = preds.poly$upr, approximating_function = "polynomial"),
  gam =  data.frame (a0 = pred_a0s, lca0 = preds.gam$mean, lwr = preds.gam$lwr, upr = preds.gam$upr, approximating_function = "gam")
)

forplot_ca0 <- do.call(rbind, preds.list)

write.csv(forplot_ca0, file = "../data/fitted_predictions_lca0_Gaussian.csv",  row.names = FALSE)

library(ggplot2)

p0 <- ggplot(data = forplot_ca0, aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = constant_data, mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  stat_function(fun = c_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0  

p0b <- ggplot(data = subset(forplot_ca0, a0 < 1), aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = subset(constant_data, a0 < 1), mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  stat_function(fun = c_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  scale_x_continuous(expression(a[0]), limits = c(0, 1)) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0b 

ggsave(p0, filename = "../figures/estimates_log_ca0_Gaussian.pdf", dpi = 300)

### Sensitivity to a_0: parameter estimates
mean_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:2, "mean"])        
  )
) 
mean_pars$a0 <- a0s[-1]
round(mean_pars, 3)
#
sd_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:2, "sd"])        
  )
) 
sd_pars$a0 <- a0s[-1]
round(sd_pars, 4)

library(ggplot2)

ggplot(data = mean_pars, aes(x = a0, y = mu)) +
  geom_line() +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(Mean(mu)), limits = c(-105, -95)) +
  geom_hline(yintercept = true.mu, linetype = "dashed") +
  theme_bw(base_size = 20) +
  NULL

ggplot(data = mean_pars, aes(x = a0, y = sigma_sq)) +
  geom_line() +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(Mean(sigma^2)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL

ggplot(data = sd_pars, aes(x = a0, y = mu)) +
  geom_line()+
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(SD(mu)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL