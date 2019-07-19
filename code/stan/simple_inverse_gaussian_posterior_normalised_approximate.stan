functions{
  int which_min(real [] y ){
    int ans = sort_indices_asc(y)[1];
    return(ans);
  }
  real approximate_ca0(real x, real[] x_pred, real[] y_pred){
    int K = size(x_pred);
    real deltas [K];
    real ans;
    int i;
    if(size(y_pred) != K) reject("x_pred and y_pred aren't of the same size");
    for(k in 1:K) deltas[k] = fabs(x_pred[k] - x);
    i = which_min(deltas);
    if(i != 1){
      real x1 = x_pred[i];
      real x2 = x_pred[i + 1];
      real y1 = y_pred[i];
      real y2 = y_pred[i + 1];
      ans = y1 + (y2-y1) * (x-x1)/(x2-x1);
    }else{
      ans = y_pred[i];
    }
    return(ans);
  }
  real inverse_gaussian_lpdf(real x, real mu, real lambda){
    real ldens = 0.5 * ( log(lambda) - log(2*pi()) - 3*log(x) ) - (lambda *(x-mu)^2)/(2* mu^2 * x); 
    return(ldens);
  }
}
data{
  int<lower=0> N0;
  real<lower=0> y0[N0];
  real<lower=0> alpha_l;
  real<lower=0> beta_l;
  real<lower=0> alpha_m;
  real<lower=0> beta_m;
  int<lower=0> N;
  real<lower=0> y[N];
  real<lower=0> delta;
  real<lower=0> nu;
  /* Approximation stuff*/
  int<lower=0> K;
  real pred_grid_x[K];
  real pred_grid_y[K];
}
parameters{
  real<lower=0> mu;
  real<lower=0> lambda;
  real<lower=0, upper=1> a_0;
}
model{
  /*Power prior*/
  target += -approximate_ca0(a_0, pred_grid_x, pred_grid_y);
  for(i in 1:N0) target += a_0 * inverse_gaussian_lpdf(y0[i] | mu, lambda);
  target += gamma_lpdf(mu | alpha_m, beta_m);
  target += gamma_lpdf(lambda | alpha_l, beta_l);
  target += beta_lpdf(a_0 | delta, nu);
  /* Likelihood */
  for(i in 1:N) target += inverse_gaussian_lpdf(y[i] | mu, lambda);
}
