data {
  int<lower=0> N0;
  int<lower=0> P;
  matrix[N0, P] X0;
  int<lower=0, upper=1> y0[N0];
  int<lower=0> N;
  matrix[N, P] X;
  int<lower=0, upper=1> y[N];
  real<lower=0> nu;
  real<lower=0> delta;
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0, upper=1> a_0;
}
model {
  /* Power prior */
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta | 0, 1);
  target += a_0 * bernoulli_logit_lpmf(y0 | alpha + X0*beta);
  target += beta_lpdf(a_0 | nu, delta);
  /* Likelihood */
  target += bernoulli_logit_lpmf(y | alpha + X*beta);
}
