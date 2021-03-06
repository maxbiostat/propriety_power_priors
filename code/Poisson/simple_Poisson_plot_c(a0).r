source("data_Poisson.r")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

po.data <- list(
  N0 = N_0,
  y0 = y_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = NA
)
###
get_l_a0_poisson <- function(y0, n0, alpha0, beta0, a_0){
  logPprime <- sum(lfactorial(y0))
  S <- sum(y0)
  ans <- -a_0 * logPprime + lgamma(a_0 * S + alpha0) - (a_0 * S + alpha0)* log(a_0 *n0 + beta0) + (alpha0*log(beta0) - lgamma(alpha0))
  return(ans)
}
############
l_a0 <- function(x) {
  get_l_a0_poisson(
    y0 = po.data$y0,
    n0 = po.data$N0,
    alpha0 = po.data$alpha0,
    beta0 = po.data$beta0,
    a_0 = x
  )
}
l_a0 <- Vectorize(l_a0)
###
l_a0_p <- function(x) numDeriv::grad(l_a0, x)
l_a0_p <- Vectorize(l_a0_p)
l_a0_pp <- function(x) numDeriv::hessian(l_a0, x)
l_a0_pp <- Vectorize(l_a0_pp)

analytical_l_p <- function(x){
  logPprime <- sum(lfactorial(y_0))
  S <- sum(y_0)
  tt <- alpha_0 + S*x
  gg <- beta_0 + N_0 * x
  result <- -logPprime -N_0 *tt/gg - S*log(gg) + S * digamma(tt)
  return(result)
}
####################
J <- 20
maxA <- 1

adaptive.ca0.estimates <- read.csv(paste("../data/constant_data/Poisson_logCA0_adaptive_J=", J, ".csv", sep = "")) 


curve(l_a0_p, 0, maxA, lwd = 2, lty = 2, col = 2,  ylab = "", xlab = expression(a[0]), ylim = c(-200, 200))
abline(h = 0 , lwd = 2 , lty = 3)
curve(l_a0, 0, maxA, lwd = 2, add = TRUE)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$lc_a0)
curve(l_a0_pp, 0, maxA, lwd = 2, col = 3, add = TRUE)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$second_deriv_lc, col = 3)
legend(x = "bottomright", legend = c("l(a_0)", "l'(a_0)", "l''(a_0)"), col = 1:3, lwd = 2, bty = 'n')
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$deriv_lc, col = 2)

######################################
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

fit.gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = adaptive.ca0.estimates)
deriv.gam.fit <- mgcv::gam(deriv_lc ~ s(a0, k = J + 1), data = adaptive.ca0.estimates)

K <- 20000
pred_a0s <- seq(0, maxA, length.out = K)

preds.gam  <- get_band(fit.gam, xpred = pred_a0s)
deriv.preds <- midpoint_integrate_la0(pred_a0s, deriv.gam.fit)

true.lca0s <- l_a0(pred_a0s)

gam.preds.list <- list(
  gam = data.frame (a0 = pred_a0s, lca0 = preds.gam$mean, dev = preds.gam$mean - true.lca0s, lwr = preds.gam$lwr,
                     upr = preds.gam$upr, approximating_function = "gam", design = "adaptive"),
  derivative = data.frame (a0 = pred_a0s, lca0 = deriv.preds, dev = deriv.preds - true.lca0s, lwr = NA, upr = NA,
                           approximating_function = "gam-derivative", design = "adaptive_numerical")
)

## RMSE
lapply(gam.preds.list, function(pred) sqrt(mean( (pred$lca0 - true.lca0s)^2 )) )
## MAD
lapply(gam.preds.list, function(pred) mean( abs(pred$lca0 - true.lca0s) )) 

## RMSE, a0 < 1
lapply(gam.preds.list, function(pred) sqrt(mean( (pred[pred$a0 < 1, ]$lca0 - true.lca0s[pred_a0s < 1])^2 )) )
## MAD, a0 < 1
lapply(gam.preds.list, function(pred) mean( abs(pred[pred$a0 < 1, ]$lca0 - true.lca0s[pred_a0s < 1]) )) 

### Exporting

forplot_ca0 <- do.call(rbind, gam.preds.list)

write.csv(forplot_ca0, file = paste("../data/constant_data/fitted_predictions_lca0_Poisson_J=", J, ".csv", sep = ""),
          row.names = FALSE)

library(ggplot2)

dev_plot <- ggplot(data = forplot_ca0, aes(x = a0, y = dev, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  scale_x_continuous(expression(a[0])) +
  theme_bw(base_size = 16)
dev_plot

p0 <- ggplot(data = forplot_ca0, aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = adaptive.ca0.estimates,
             mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  stat_function(fun = l_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0  

ggsave(p0, filename = paste("../figures/estimates_log_ca0_Poisson_J=", J, ".pdf", sep = ""), dpi = 300)
