library(npowerPrioR)
scenario <- 1
source(paste("data_regression_NIG_scenario_", scenario, ".r", sep = ""))
summary(lm(y_0 ~ -1 + X_0))

lm.data <- list(
  N0 = N_0,
  P = P,
  X0 = X_0,
  y0 = y_0,
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  a_0 = NULL
)

invlambda0 <- solve(lm.data$lambda_0)
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
l_a0 <- function(x) {
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
l_a0 <- Vectorize(l_a0)
l_a0_p <- function(x) numDeriv::grad(l_a0, x)
l_a0_p <- Vectorize(l_a0_p)
l_a0_pp <- function(x) numDeriv::hessian(l_a0, x)
l_a0_pp <- Vectorize(l_a0_pp)

#########
J <- 20

adaptive.ca0.estimates <- read.csv(paste("../data/constant_data/RegressionNIG_logCA0_adaptive_J=", J,
                                         "_scenario_", scenario, ".csv", sep = ""))
maxA <- max(adaptive.ca0.estimates$a0)
ylims <- range(c(adaptive.ca0.estimates$lc_a0, adaptive.ca0.estimates$deriv_lc))

curve(l_a0_p, 1E-4, maxA, lwd = 3, lty = 2, col = 2, 
      ylab = "", xlab = expression(a[0]), cex.lab = 3, cex.axis = 2.5, ylim = ylims)
abline(h = 0, lwd = 2 , lty = 3)
curve(l_a0, 0, maxA, lwd = 3, add = TRUE)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$lc_a0, cex = 1.5, pch = 16)
curve(l_a0_pp, 0, maxA, lwd = 3, col = 3, add = TRUE)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$second_deriv_lc, cex = 1.5, col = 3, pch = 16)
points(adaptive.ca0.estimates$a0, adaptive.ca0.estimates$deriv_lc,
       pch = 16, cex = 1.5, col = 2)
legend(x = "topright", legend = c("l(a_0)", "l'(a_0)", "l''(a_0)"),
       col = 1:3, lwd = 2, bty = 'n', cex = 1.5)


adaptive.ca0.estimates$true_la0 <- l_a0(adaptive.ca0.estimates$a0)
adaptive.ca0.estimates$norm_error <- (adaptive.ca0.estimates$lc_a0- adaptive.ca0.estimates$true_la0)/adaptive.ca0.estimates$true_la0

adaptive.ca0.estimates
mean(abs(adaptive.ca0.estimates$norm_error[-1]))


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

fit.gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = adaptive.ca0.estimates)
deriv.gam.fit <- mgcv::gam(deriv_lc ~ s(a0, k = J + 1), data = adaptive.ca0.estimates)

K <- 20000
pred_a0s <- seq(0, maxA, length.out = K)

preds.gam  <- get_band(fit.gam, xpred = pred_a0s)
deriv.preds <- midpoint_integrate_la0(pred_a0s, deriv.gam.fit)

true.lca0s <- l_a0(pred_a0s)

gam.preds.list <- list(
  gam =  data.frame (a0 = pred_a0s, lca0 = preds.gam$mean, dev = preds.gam$mean - true.lca0s, lwr = preds.gam$lwr,
                     upr = preds.gam$upr, approximating_function = "gam", design = "adaptive"),
  derivative = data.frame (a0 = pred_a0s, lca0 = deriv.preds, dev = deriv.preds - true.lca0s, lwr = NA, upr = NA,
                           approximating_function = "gam-derivative", design = "adaptive_numerical")
)

##########################

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

write.csv(forplot_ca0,
          file = paste("../data/constant_data/fitted_predictions_lca0_RegressionNIG_scenario_", scenario, ".csv", sep = ""),  row.names = FALSE)

library(ggplot2)

p0 <- ggplot(data = forplot_ca0, aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = adaptive.ca0.estimates, mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  stat_function(fun = l_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom",
        legend.justification = "centre",
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0))
p0  

p1 <- ggplot(data = subset(forplot_ca0, approximating_function == "gam"),
             aes(x = a0, y = lca0, colour = approximating_function, fill = approximating_function)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = adaptive.ca0.estimates, mapping = aes(x = a0, y = lc_a0), colour = "black", alpha = .4, size = 5, inherit.aes = FALSE) + 
  stat_function(fun = l_a0, linetype = "dashed", size = 1.2, colour = "black") + 
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0])))) +
  theme_bw(base_size = 20) + 
  theme(legend.position = "bottom",
        legend.justification = "centre",
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0))
p1 

ggsave(p0,
  filename = paste("../figures/estimates_log_ca0_RegressionNIG_scenario_", scenario, ".pdf", sep = ""), dpi = 300)