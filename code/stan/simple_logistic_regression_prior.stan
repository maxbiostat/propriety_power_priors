data {
  int<lower=0> N0;
  int<lower=0> P;
  matrix[N0, P] X0;
  int<lower=0, upper=1> y0[N0];
  real<lower=0> a_0;
}
parameters {
  real alpha;
  vector[P] beta;
}
model {
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta | 0, 1);
  target += a_0 * bernoulli_logit_lpmf(y0 | alpha + X0*beta);
}
