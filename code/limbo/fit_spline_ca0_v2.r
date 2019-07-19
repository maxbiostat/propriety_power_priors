library(splines)
library(rstan)
library(ggplot2)
library(MonoPoly)

extra_width <- 0.1
x_pred <- seq(from = -1 - extra_width , to = 1 + extra_width, length.out = 1000)

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
    hawkins$x,
    knots = internal_knots,
    Boundary.knots = boundary_knots
  )
  
  bs_pred_mat <- splines2::cSpline(
    x_pred,
    knots = internal_knots,
    Boundary.knots = boundary_knots
  )
}else{
  bs_data_mat <- bs(
    hawkins$x,
    knots = internal_knots,
    Boundary.knots = boundary_knots
  )
  
  bs_pred_mat <- bs(
    x_pred,
    knots = internal_knots,
    Boundary.knots = boundary_knots
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
  positive_ordered [n_knots] theta;
  real <lower = 0> y_sd;
  real intercept;
}

model {
  theta ~ normal(0, 15);
  y ~ normal(intercept + x_data_mat * theta, y_sd);
  y_sd ~ normal(0, 4);
  intercept ~ normal(0, 10);
}

generated quantities {
  vector [n_pred_points] post_pred_y = intercept + x_pred_mat * theta;
}
'
# 
stan_data <- list(
  n_knots = ncol(bs_data_mat),
  n_pred_points = length(x_pred),
  x_pred_mat = bs_pred_mat,
  n_data_points = nrow(hawkins),
  x_data_mat = bs_data_mat,
  y = hawkins$y
)

model_prefit <- stan_model(model_code = spline_code)
model_fit <- sampling(
  model_prefit,
  data = stan_data,
  cores = 4
)

post_pred_samples <- rstan::extract(model_fit, "post_pred_y")[[1]]
post_pred_quantiles <- matrixStats::colQuantiles(post_pred_samples, probs = c(0.01, 0.5, 0.99))

plot_df <- data.frame(
  x = x_pred,
  lower = post_pred_quantiles[, 1],
  median = post_pred_quantiles[, 2],
  upper = post_pred_quantiles[, 3]
)

ggplot(plot_df, aes(x = x)) +
  geom_point(data = hawkins, inherit.aes = FALSE, aes(x = x, y = y)) + 
  geom_line(aes(y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25)