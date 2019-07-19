source("power_priors_aux.r")
source("Bernoulli_data.r")

library(rstan)
rstan_options(auto_write = TRUE)

options(mc.cores = 4)
rstan:::rstudio_stanc("stan/simple_bernoulli_prior.stan")

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_bernoulli_prior.stan")

bb.data <- list(
  N0 = N_0,
  y0 = y_0,
  c = cc,
  d = dd,
  a_0 = .4
)

get_c_a0_bernoulli <- function(y0, n0, cc, dd, a_0){
  ans <- lbeta(a_0 * y0 + cc, a_0 *(n0 -y0) + dd)
  return(ans)
}

get_c_a0_bernoulli(
  y0 = bb.data$y0,
  n0 = bb.data$N0,
  cc = bb.data$c,
  dd = bb.data$d,
  a_0 = bb.data$a_0
)

anaughts <- seq(0, 2, length.out = 100)
c_a0 <- function(x) {
  get_c_a0_bernoulli(
    y0 = bb.data$y0,
    n0 = bb.data$N0,
    cc = bb.data$c,
    dd = bb.data$d,
    a_0 = x
  )
} 

#####################
get_estimates_a0 <- function(x){
  bb.data <- list(
    N0 = N_0,
    y0 = y_0,
    c = cc,
    d = dd,
    a_0 = x
  )
  res <- list(chain = sampling(compiled.model.prior, data = bb.data, refresh = 0))
  res$log_mal <- bridgesampling::bridge_sampler(res$chain, silent = TRUE)            
  return(res)
}
J <- 10
maxA <- 5
# est_a0s <- seq(0.05, maxA, length.out = J)
# est_a0s <- pracma::logseq(0.01, maxA, n = J)
est_a0s <- c(seq(0.05, 1, length.out = J-3), seq(1.2, maxA, length.out = 3))
system.time(
  results <- lapply(est_a0s, get_estimates_a0)  
)

######################################
constant_data <- data.frame(
  a0 = est_a0s,
  lc_a0 = unlist(lapply(results, function(r) r$log_mal$logml))
)

write.csv(constant_data, file = "../data/Bernoulli_logCA0.csv", row.names = FALSE)

plot(constant_data)

degree <- 2

fit.poly <- lm(lc_a0 ~ poly(a0, degree), data = constant_data)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0), data = constant_data)


K <- 50
pred_a0s <- seq(0, maxA, length.out = K)

preds.poly <- get_band(fit.poly, xpred = pred_a0s)
preds.gam  <- get_band(fit.gam, xpred = pred_a0s)

preds.list <- list(
  polynomial = data.frame (a0 = pred_a0s, lca0 = preds.poly$mean,
                           lwr = preds.poly$lwr, upr = preds.poly$upr, approximating_function = "polynomial"),
  gam =  data.frame (a0 = pred_a0s, lca0 = preds.gam$mean,
                     lwr = preds.gam$lwr, upr = preds.gam$upr, approximating_function = "gam")
)

forplot_ca0 <- do.call(rbind, preds.list)

write.csv(forplot_ca0, file = "../data/fitted_predictions_lca0_Bernoulli.csv",  row.names = FALSE)

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

ggsave(p0, filename = "../figures/estimates_log_ca0_Bernoulli.pdf", dpi = 300)

p0b <- ggplot(data = subset(forplot_ca0, a0 < 1), aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = subset(constant_data, a0 < 1), mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  stat_function(fun = c_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  scale_x_continuous(expression(a[0]), limits = c(0, 1)) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0b 
######################################
### Sensitivity to a_0: parameter estimates
mean_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:2, "mean"])        
  )
) 
mean_pars$a0 <- est_a0s
round(mean_pars, 3)
#
sd_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:2, "sd"])        
  )
) 
sd_pars$a0 <- est_a0s
round(sd_pars, 4)

library(ggplot2)

ggplot(data = mean_pars, aes(x = a0, y = theta)) +
  geom_line() +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(Mean(theta)), expand = c(0, 0), limits = c(.075 , .15)) +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = .1, linetype = "dashed") +
  NULL

ggplot(data = sd_pars, aes(x = a0, y = theta)) +
  geom_line()+
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(SD(theta)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL