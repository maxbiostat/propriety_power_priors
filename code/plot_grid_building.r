J  <- 10
adaptive.ca0.estimates <- read.csv(paste("../data/constant_data/Gaussian_logCA0_adaptive",
                                         "_J=", J, ".csv", sep = ""))

adaptive.ca0.estimates$order <- c(0, 1:(nrow(adaptive.ca0.estimates)-2), 1)

fit.gam.adaptive <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = adaptive.ca0.estimates)

K <- 20000
maxA <- max(adaptive.ca0.estimates$a0)
pred_a0s <- seq(0, maxA, length.out = K)

preds.gam.adaptive  <- get_band(fit.gam.adaptive, xpred = pred_a0s)

adaptive.preds <-  data.frame (a0 = pred_a0s, lca0 = preds.gam.adaptive$mean,
                               lwr = preds.gam.adaptive$lwr,
                        upr = preds.gam.adaptive$upr)

library(ggplot2)

selected <- subset(adaptive.ca0.estimates, a0 > 0 & order < 6) #
selected$point_name <- paste0("j=", selected$order)
selected$lca0 <- selected$lc_a0

p0 <- ggplot(data = adaptive.preds, aes(x = a0, y = lca0)) +
  geom_line() +
  geom_ribbon(aes(min = lwr, max = upr), alpha = .4) +
  geom_point(data = selected ,
             mapping = aes(x = a0, y = lc_a0), alpha = .5, size = 4,
             alpha = .8, inherit.aes = FALSE) +
  geom_text(data = selected, mapping = aes(label = point_name), hjust = 0, vjust = -2) +
  scale_x_continuous(expression(a[0])) +
  scale_y_continuous(expression(log(c(a[0]))), limits = c(-6, 690)) +
  theme_bw(base_size = 16)
p0  

ggsave(p0, filename = "../figures/grid_building_schema.pdf", dpi = 500)

