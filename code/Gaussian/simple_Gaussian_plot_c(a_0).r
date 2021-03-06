library(npowerPrioR)
source("data_Gaussian.r")

gs.data <- list(
  N0 = N_0,
  y0 = y_0,
  mu0 = mu_0,
  kappa0 = kappa_0,
  alpha0 = alpha_0,
  beta0 = beta_0,
  a_0 = 1
)

###
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
    y0 = gs.data$y0,
    n0 = gs.data$N0,
    alpha0 = gs.data$alpha0,
    beta0 = gs.data$beta0,
    m0 = gs.data$mu0,
    k0 = gs.data$kappa0,
    a_0 = x
  )
}
l_a0 <- Vectorize(l_a0)
###
l_a0_p <- function(x) numDeriv::grad(l_a0, x)
l_a0_p <- Vectorize(l_a0_p)
l_a0_pp <- function(x) numDeriv::hessian(l_a0, x)
l_a0_pp <- Vectorize(l_a0_pp)

####################
J <- 20
maxA <- 10

## Importing 
uniform.ca0.estimates <- read.csv(paste("../data/constant_data/Gaussian_logCA0_uniform", "_J=", J, ".csv", sep = ""))
adaptive.ca0.estimates <- read.csv(paste("../data/constant_data/Gaussian_logCA0_adaptive", "_J=", J, ".csv", sep = ""))
adaptive.ca0.estimates.derivOnly <- read.csv(paste("../data/constant_data/Gaussian_logCA0_adaptive_derivOnly", "_J=", J, ".csv", sep = ""))

####
## Plotting the curves l(a_0) and derivatives

curve(l_a0_p, 0, 10, lwd = 3, lty = 2, col = 2, 
      # ylim = c(-30, 60),
      ylab = "", xlab = expression(a[0]), cex.lab = 3, cex.axis = 2.5)
abline(h = 0, lwd = 2 , lty = 3)
curve(l_a0, 0, 10, lwd = 3, add = TRUE)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$lc_a0, cex = 1.5, pch = 16)
curve(l_a0_pp, 0, 10, lwd = 3, col = 3, add = TRUE)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$second_deriv_lc, cex = 1.5, col = 3, pch = 16)
legend(x = "bottomright", legend = c("l(a_0)", "l'(a_0)", "l''(a_0)"),
       col = 1:3, lwd = 2, bty = 'n', cex = 1.5)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$deriv_lc,
       pch = 16, cex = 1.5, col = 2)


#######

fit.gam.uniform <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = uniform.ca0.estimates)
fit.gam.adaptive <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = adaptive.ca0.estimates)
deriv.gam.fit <- mgcv::gam(deriv_lc ~ s(a0, k = J + 1), data = adaptive.ca0.estimates.derivOnly)

########################
K <- 20000
pred_a0s <- seq(0, maxA, length.out = K)

preds.gam.uniform  <- get_band(fit.gam.uniform, xpred = pred_a0s)
preds.gam.adaptive  <- get_band(fit.gam.adaptive, xpred = pred_a0s)
deriv.preds <- cumsum(predict(deriv.gam.fit, newdata = data.frame(a0 = pred_a0s))) * diff(pred_a0s)[1]

true.lca0s <- l_a0(pred_a0s)
gam.preds.list <- list(
  uniform =  data.frame (a0 = pred_a0s, lca0 = preds.gam.uniform$mean, lwr = preds.gam.uniform$lwr,
                         upr = preds.gam.uniform$upr, approximating_function = "gam", design = "uniform"),
  adaptive =  data.frame (a0 = pred_a0s, lca0 = preds.gam.adaptive$mean, lwr = preds.gam.adaptive$lwr,
                          upr = preds.gam.adaptive$upr, approximating_function = "gam", design = "adaptive"),
  derivative = data.frame (a0 = pred_a0s, lca0 = deriv.preds, lwr = NA, upr = NA,
                           approximating_function = "numerical", design = "adaptive_numerical")
)


l_a0(.15)
predict(fit.gam.uniform, newdata = data.frame(a0 = .15))
predict(fit.gam.adaptive, newdata = data.frame(a0 = .15))

l_a0(.75)
predict(fit.gam.uniform, newdata = data.frame(a0 = .75))
predict(fit.gam.adaptive, newdata = data.frame(a0 = .75))

l_a0(5)
predict(fit.gam.uniform, newdata = data.frame(a0 = 5))
predict(fit.gam.adaptive, newdata = data.frame(a0 = 5)) 

## RMSE
lapply(gam.preds.list, function(pred) sqrt(mean( (pred$lca0 - true.lca0s)^2 )) )
## MAD
lapply(gam.preds.list, function(pred) mean( abs(pred$lca0 - true.lca0s) )) 

## RMSE, a0 < 1
lapply(gam.preds.list, function(pred) sqrt(mean( (pred[pred$a0 < 1, ]$lca0 - true.lca0s[pred_a0s < 1])^2 )) )
## MAD, a0 < 1
lapply(gam.preds.list, function(pred) mean( abs(pred[pred$a0 < 1, ]$lca0 - true.lca0s[pred_a0s < 1]) )) 

forplot_ca0 <- do.call(rbind, gam.preds.list)

write.csv(forplot_ca0,
          file = paste("../data/constant_data/fitted_predictions_lca0_Gaussian", "_J=", J, ".csv", sep = ""),
          row.names = FALSE)

library(ggplot2)

p0 <- ggplot(data = forplot_ca0, aes(x = a0, y = lca0, colour = design, fill = design)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  stat_function(fun = l_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  geom_point(data = adaptive.ca0.estimates,
             mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0  

p0b <- ggplot(data = subset(forplot_ca0, a0 < 1 & design != "uniform"), aes(x = a0, y = lca0, colour = design, fill = design)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  stat_function(fun = l_a0, linetype = "dashed", size = 1.2, colour = "black") +
  geom_point(data = subset(adaptive.ca0.estimates, a0 < 1),
             mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  scale_x_continuous(expression(a[0]), limits = c(0, 1)) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 16)
p0b

plot(pred_a0s, (true.lca0s-preds.gam.adaptive$mean))
plot(pred_a0s, (true.lca0s-deriv.preds))

ggsave(p0, filename = paste("../figures/estimates_log_ca0_Gaussian_J=", J, ".pdf", sep = ""), dpi = 300)

ggsave(p0b, filename = paste("../figures/estimates_log_ca0_Gaussian_restricted_J=", J, ".pdf", sep = ""), dpi = 300)