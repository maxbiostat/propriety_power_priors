source("logistic_regression_data.r")
source("power_priors_aux.r")
summary(glm(y_0 ~  X_0, family = "binomial"))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"

prior.model <- stan_model("stan/simple_logistic_regression_prior.stan")

# lgr.data <- list(
#   N0 = N_0,
#   X0 = X_0,
#   y0 = y_0,
#   a_0 = NULL
# )

get_c_a0_brute_force <- function(x){
  lgr.data <- list(
    N0 = N_0,
    X0 = X_0,
    y0 = y_0,
    a_0 = x
  )
  prior.lgr <- sampling(prior.model, data = lgr.data, refresh = 0)
  res <- list(
    chain = prior.lgr,
    logml = bridgesampling::bridge_sampler(prior.lgr, silent = TRUE)$logml
  )
  return(res)
}
###
J <- 20
maxA <- 10
# est_a0s <- seq(0.1, maxA, length.out = J)
# est_a0s <- pracma::logseq(0.01, maxA, n = J)
est_a0s <- c(seq(0.01, 1, length.out = J-3), seq(1.2, maxA, length.out = 3))
system.time(
 results <- lapply(est_a0s, get_c_a0_brute_force)  
)
mls <- unlist(lapply(results, function(x) x$logml))

constant_data <- data.frame(a0 = c(0, est_a0s), lc_a0 = c(0, mls))

write.csv(constant_data, "../data/RegressionLogistic_logCA0.csv", row.names = FALSE)

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
  polynomial = data.frame (a0 = pred_a0s , lca0 = preds.poly$mean, lwr = preds.poly$lwr,
                           upr = preds.poly$upr, approximating_function = "polynomial"),
  gam =  data.frame (a0 = pred_a0s, lca0 = preds.gam$mean, lwr = preds.gam$lwr, upr = preds.gam$upr, approximating_function = "gam")
)

forplot_ca0 <- do.call(rbind, preds.list)

write.csv(forplot_ca0, file = "../data/fitted_predictions_lca0_RegressionLogistic.csv",  row.names = FALSE)

library(ggplot2)

p0 <- ggplot(data = forplot_ca0, aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = constant_data, mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0  

ggsave(p0, filename = "../figures/estimates_log_ca0_RegressionLogistic.pdf", dpi = 300)

p0b <- ggplot(data = subset(forplot_ca0, a0 < 1),
              aes(x = a0, y = lca0,
                  colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = subset(constant_data, a0 < 1),
             mapping = aes(x = a0, y = lc_a0),
             colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  scale_x_continuous(expression(a[0]), limits = c(0, 1)) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0b 

### Sensitivity to a_0: parameter estimates
mean_betas <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[2:5, "mean"])        
  )
) 
mean_betas$a0 <- est_a0s
round(mean_betas, 3)
#
sd_betas <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:5, "sd"])        
  )
) 
sd_betas$a0 <- est_a0s
round(sd_betas, 3)
