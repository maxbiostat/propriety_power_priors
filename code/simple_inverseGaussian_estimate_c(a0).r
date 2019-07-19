source("power_priors_aux.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_inverse_gaussian_prior.stan")

source("inverse_Gaussian_data.r")
summary(y0)

ig.data <- list(
  N0 = N0,
  y0 = y0,
  alpha_l = a.l,
  beta_l = b.l,
  alpha_m = a.m,
  beta_m = b.m,
  a_0 = NULL
)

get_c_a0_brute_force <- function(x){
  ig.data <- list(
    N0 = N0,
    y0 = y0,
    alpha_l = a.l,
    beta_l = b.l,
    alpha_m = a.m,
    beta_m = b.m,
    a_0 = x
  )
  prior.invgauss <- sampling(compiled.model.prior, data = ig.data, refresh = 0)
  res <- list(
    chain =  prior.invgauss,
    logml = bridgesampling::bridge_sampler(prior.invgauss, silent = TRUE)$logml
  )
  return(res)
}

###
get_c_a0_invgauss <- function(y0, n0, alpha.l, beta.l, alpha.m, beta.m, a_0){
  invgaus_kernel_forint_2 <- function(u, a, b, c, d, s, r){
    (a/(c * u^2) - b/(c* u) + 1)^-d * u^(s-1) * exp(-r * u)
  }
  invgaus_kernel_forint_2 <- Vectorize(invgaus_kernel_forint_2)
  ###########
  logP <- sum(log(y0))
  S <- sum(y0)
  Sprime <- sum(1/y0)
  C <- a_0*Sprime/2 + beta.l
  D <- (a_0*n0)/2 + alpha.l
  lconst <- -0.5 *n0*a_0*log(2*pi) - (3*a_0/2)*logP + lgamma(D) + alpha.m * log(beta.m) +  alpha.l * log(beta.l) - lgamma(alpha.l) - lgamma(alpha.m)
  kf <- function(x, 
                 a = a_0*S/2,
                 b = a_0*n0,
                 c = C,
                 d = D,
                 s = alpha.m,
                 r = beta.m){
    invgaus_kernel_forint_2(u = x, a = a, b = b, c = c, d = d, s = s, r = r)
  }
  Int <- tryCatch(integrate(kf, 0, Inf), error = function(e) return(list(value = NA)))
  ###
  ans <- lconst  -D * log(C) +  log(Int$value)
  return(ans)
}
############
c_a0 <- function(x) {
  get_c_a0_invgauss(
    y0 = ig.data$y0,
    n0 = ig.data$N0,
    alpha.l = ig.data$alpha_l,
    beta.l = ig.data$beta_l,
    alpha.m = ig.data$alpha_m,
    beta.m = ig.data$beta_m,
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
write.csv(constant_data, row.names = FALSE, "../data/InvGaussian_logCA0.csv")

#######
K <- 10000
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

write.csv(forplot_ca0, file = "../data/fitted_predictions_lca0_invgaussian.csv",  row.names = FALSE)

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

ggsave(p0, filename = "../figures/estimates_log_ca0_invGaussian.pdf", dpi = 300)

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
  scale_y_continuous(expression(Mean(mu))) +
  geom_hline(yintercept = true.mu, linetype = "dashed") +
  theme_bw(base_size = 20) +
  NULL

ggplot(data = mean_pars, aes(x = a0, y = lambda)) +
  geom_line() +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(Mean(lambda)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL

ggplot(data = sd_pars, aes(x = a0, y = mu)) +
  geom_line()+
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(SD(mu)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL

ggplot(data = sd_pars, aes(x = a0, y = lambda)) +
  geom_line()+
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(SD(lambda)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL
