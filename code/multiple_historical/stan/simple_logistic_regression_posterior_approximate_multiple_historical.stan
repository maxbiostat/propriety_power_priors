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
}
data {
  int<lower=0> N0max;
  int<lower=0> D;
  int<lower=0> N0[D];
  int<lower=0> P;
  real X0[N0max, P, D];
  int<lower=0, upper=1> y0[N0max, D];
  int<lower=0> N;
  matrix[N, P] X;
  int<lower=0, upper=1> y[N];
  real<lower=0> eta[D];
  real<lower=0> nu[D];
  /* Approximation stuff*/
  int<lower=0> K;
  real pred_grid_x[K, D];
  real pred_grid_y[K, D];
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0, upper=1> a_0[D];
}
model {
  /* Power prior */
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta | 0, 1);
  for(d in 1:D){
  target += -approximate_ca0(a_0[d], pred_grid_x[, d], pred_grid_y[, d]);
  target += a_0[d] * bernoulli_logit_lpmf( y0[1:N0[d], d] | alpha + to_matrix(X0[1:N0[d], , d])*beta);
  target += beta_lpdf(a_0[d] | eta[d], nu[d]);
  }
  /* Likelihood */
  target += bernoulli_logit_lpmf(y | alpha + X*beta);
}
