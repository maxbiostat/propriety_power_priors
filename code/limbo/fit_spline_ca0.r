library(splines)
library(rstan)
library(ggplot2)
library(MonoPoly)

ca0.data <- read.csv("Gaussian_logCA0.csv")

x_pred <- seq(from = 0, to = 5, length.out = 50)

n_knots <- 10
boundary_knots <- c(min(x_pred), max(x_pred))
internal_knots <- seq(
  from = boundary_knots[1],
  to = boundary_knots[2],
  length.out = n_knots
)[2 : (n_knots - 1)]


convex <- TRUE

if(convex){
  bs_data_mat <- splines2::cSpline(
    x = ca0.data$a0,
    df = 5,
    knots = internal_knots,
    Boundary.knots = boundary_knots
  )
  bs_pred_mat <- splines2::cSpline(
    x_pred,
    df = 5,
    knots = internal_knots,
    Boundary.knots = boundary_knots
  )
}else{
  bs_data_mat <- bs(
    x = ca0.data$a0,
    df = 5,
    knots = internal_knots,
    Boundary.knots = boundary_knots,
    intercept = FALSE
  )
  bs_pred_mat <- bs(
    x_pred,
    df = 5,
    knots = internal_knots,
    Boundary.knots = boundary_knots,
    intercept = FALSE
  )
}

#
spline_code <- '
data {
  int n_knots;
  
  // the prediction spline
  int n_pred_points;
  matrix [n_pred_points, n_knots] x_pred_mat;
  
  // the data spline
  int n_data_points;
  matrix [n_data_points, n_knots] x_data_mat;
  
  vector [n_data_points] y;
}

parameters {
  vector[n_knots] theta;
  real <lower = 0> y_sd;
}

model {
  theta ~ normal(0, 15);
  y ~ normal(x_data_mat * theta, y_sd);
  y_sd ~ normal(0, 4);
}

generated quantities {
  vector [n_pred_points] post_pred_y = x_pred_mat * theta;
}
'
# 
stan_data <- list(
  n_knots = ncol(bs_data_mat),
  n_pred_points = length(x_pred),
  x_pred_mat = bs_pred_mat,
  n_data_points = nrow(ca0.data),
  x_data_mat = bs_data_mat,
  y = ca0.data$lc_a0
)

model_prefit <- stan_model(model_code = spline_code)
model_fit <- sampling(
  model_prefit,
  data = stan_data,
  cores = 4
)

model_fit

post_pred_samples <- rstan::extract(model_fit, "post_pred_y")[[1]]
post_pred_quantiles <- matrixStats::colQuantiles(post_pred_samples, probs = c(0.01, 0.5, 0.99))

library(brms)
m0 <- brm(lc_a0 ~s(a0), data = ca0.data, control = list(adapt_delta = .99))

plot(m0)

post_pred_quantiles <- predict(m0, newdata = data.frame(a0 = x_pred))

plot_df <- data.frame(
  x = x_pred,
  lower = post_pred_quantiles[, 3],
  median = post_pred_quantiles[, 1],
  upper = post_pred_quantiles[, 4]
)

ggplot(plot_df, aes(x = x)) +
  geom_point(data = ca0.data, inherit.aes = FALSE, aes(x = a0, y = lc_a0)) +
  geom_line(aes(y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25)
