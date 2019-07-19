
source("regression_NIG_data.r")
source("power_priors_aux.r")
summary(lm(y_0 ~ -1 + X_0))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#### Sampling from the "prior"

prior.model <- stan_model("stan/simple_linear_regression_NIG_prior.stan")

lm.data <- list(
  N0 = N_0,
  X0 = X_0,
  y0 = y_0,
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  a_0 = NULL
)

invlambda0 <- solve(lm.data$lambda_0)
get_c_a0_brute_force <- function(x){
  lm.data <- list(
    N0 = N_0,
    X0 = X_0,
    y0 = y_0,
    mu_beta = rep(0, P),
    lambda_0 = solve(vb * diag(P)),
    alpha0 = as,
    beta0 = bs,
    a_0 = x
  )
  prior.lm <- sampling(prior.model, data = lm.data, refresh = 0)
  res <- list(
    chain = prior.lm,
    logml = bridgesampling::bridge_sampler(prior.lm, silent = TRUE)$logml
  )
  return(res)
}
###
get_mal_NIG_regression <- function(y0, X0, n0, mu0, lambda0, alpha0, beta0, a_0){
  P <- ncol(X0)
  if(length(mu0) != P) stop("mu0 is not the same dimension as X")
  Xstar <- sqrt(a_0) * X0
  ystar <- sqrt(a_0) * y0
  lambda_n <- t(Xstar)%*%Xstar + lambda0   
  mu_n <- solve(lambda_n) %*% (lambda0%*%mu0 + t(Xstar)%*%ystar)
  alpha_n <-  as + (n0*a_0)/2
  beta_n <- beta0 + .5 * ( t(ystar)%*%ystar + t(mu0)%*%lambda0%*%mu0 - t(mu_n)%*%lambda_n%*%mu_n )
  det0 <- det(lambda0)
  det1 <- det(lambda_n)   
  ans <- -(a_0 * n0/2) * log(2*pi) + .5 * (log(det0) - log(det1)) + (alpha0 * log(beta0) - alpha_n * log(beta_n) ) + (lgamma(alpha_n)-lgamma(alpha0))  
  return(ans)
}
#
c_a0 <- function(x) {
  get_mal_NIG_regression(
    y0 = lm.data$y0,
    X0 = lm.data$X0,
    n0 = lm.data$N0,
    mu0 = lm.data$mu_beta,
    lambda0 = invlambda0,
    alpha0 = lm.data$alpha0,
    beta0 = lm.data$beta0,
    a_0 = x
  )
}
c_a0 <- Vectorize(c_a0)
##
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

write.csv(constant_data, "../data/RegressionNIG_logCA0.csv", row.names = FALSE)

#######
K <- 10000
pred_a0s <- seq(0, maxA, length.out = K)
degree <- 2

fit.poly <- lm(lc_a0 ~ poly(a0, degree), data = constant_data)
fit.gam <- mgcv::gam(lc_a0 ~ s(a0), data = constant_data)

predict(fit.gam, newdata = data.frame(a0 = .3))
c_a0(.3)

predict(fit.gam, newdata = data.frame(a0 = .03))
c_a0(.03)

##########################
preds.poly <- get_band(fit.poly, xpred = pred_a0s)
preds.gam  <- get_band(fit.gam, xpred = pred_a0s)

preds.list <- list(
  polynomial = data.frame (a0 = pred_a0s , lca0 = preds.poly$mean, lwr = preds.poly$lwr,
                           upr = preds.poly$upr, approximating_function = "polynomial"),
  gam =  data.frame (a0 = pred_a0s, lca0 = preds.gam$mean, lwr = preds.gam$lwr, upr = preds.gam$upr, approximating_function = "gam")
)

forplot_ca0 <- do.call(rbind, preds.list)

write.csv(forplot_ca0, file = "../data/fitted_predictions_lca0_RegressionNIG.csv",  row.names = FALSE)

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

ggsave(p0, filename = "../figures/estimates_log_ca0_RegressionNIG.pdf", dpi = 300)

p0b <- ggplot(data = subset(forplot_ca0, a0 < 1),
              aes(x = a0, y = lca0,
                  colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = subset(constant_data, a0 < 1),
             mapping = aes(x = a0, y = lc_a0),
             colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  stat_function(fun = c_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  scale_x_continuous(expression(a[0]), limits = c(0, 1)) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0b 

### Sensitivity to a_0: parameter estimates
mean_betas <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r$chain)$summary[1:4, "mean"])        
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
