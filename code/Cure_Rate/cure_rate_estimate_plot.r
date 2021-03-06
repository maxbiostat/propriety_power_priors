source("power_priors_aux.r")
source("data_cure_rate_real.r")

#### Sampling from the "prior"

cr.data <- list(
  n = N_0,
  p = ncol(X_0),
  y = Y_0_cens,
  X = X_0,
  delta = Delta_0,
  mu_0 = mu0,
  sigma_0 = s0,
  delta_0 = d0,
  tau_0 = tau0,
  sigma_beta = 10,
  a_0 = NA
)
####################
J <- 20

load("../data/sensitivity_data/cure_rate_real.RData")

constant.data <- adaptive.ca0.estimates$result

maxA <- max(constant.data$a0)

#######

midpoint_integrate_la0 <- function(x_preds, fit){
  ## takes a grid of values for the input and a fit object that gives predictions (GLM, GAM, etc) and
  ## Gives the integrated function
  n <- length(x_preds)
  delta <- (x_preds[2]-x_preds[1])
  midpoints <- 0.5*(x_preds[2:n] + x_preds[1:(n-1)])
  ans <-  cumsum(predict(fit, newdata = data.frame(a0 = midpoints))) * delta
  ans <- c(0, ans)
  return(ans)
}

fit.gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant.data)
deriv.gam.fit <- mgcv::gam(deriv_lc ~ s(a0, k = J), data = constant.data[-1, ])

K <- 20000
pred_a0s <- seq(0, maxA, length.out = K)

preds.gam  <- get_band(fit.gam, xpred = pred_a0s)
deriv.preds <- midpoint_integrate_la0(pred_a0s, deriv.gam.fit)


gam.preds.list <- list(
  gam =  data.frame (a0 = pred_a0s, lca0 = preds.gam$mean, lwr = preds.gam$lwr,
                     upr = preds.gam$upr, approximating_function = "gam", design = "adaptive"),
  derivative = data.frame (a0 = pred_a0s, lca0 = deriv.preds, lwr = NA, upr = NA,
                           approximating_function = "gam-derivative", design = "adaptive_numerical")
)


forplot_ca0 <- do.call(rbind, gam.preds.list)

write.csv(forplot_ca0, file = "../data/constant_data/fitted_predictions_lca0_CureRateSimuData.csv",  row.names = FALSE)

### Plotting 

library(ggplot2)

p0 <- ggplot(data = forplot_ca0, aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = constant.data, mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0  


ggsave(p0, filename = "../figures/estimates_log_ca0_CureRateRealData.pdf", dpi = 300)

####
## Sensitivity analysis
sens.analysis.dt <- do.call(rbind, lapply(1:length(adaptive.ca0.estimates$summaries), function(i) {
  y <- data.frame(adaptive.ca0.estimates$summaries[[i]])
  y$a_0 <- adaptive.ca0.estimates$result$a0[i]
  return(y)
})
)

p1 <- ggplot(data = sens.analysis.dt, aes(x = a_0, y = mean, fill = parameter, colour = parameter)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0), limit = c(0, 1)) +
  scale_y_continuous("", expand = c(0, 0)) +
  facet_grid(parameter~., labeller = label_parsed, scales = "free")+
  theme_bw(base_size = 16)
p1

ggsave(p1, filename = "../figures/sensitivity_CureRateRealData.pdf", dpi = 300)

p1b <- ggplot(data = subset(sens.analysis.dt, a_0  < 0.2), aes(x = a_0, y = mean, fill = parameter, colour = parameter)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous("", expand = c(0, 0)) +
  facet_grid(parameter~., labeller = label_parsed, scales = "free")+
  theme_bw(base_size = 16)
p1b

p1c <- ggplot(data = subset(sens.analysis.dt, a_0 > 0.2), aes(x = a_0, y = mean, fill = parameter, colour = parameter)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous("", expand = c(0, 0)) +
  facet_grid(parameter~., labeller = label_parsed, scales = "free")+
  theme_bw(base_size = 16)
p1c
