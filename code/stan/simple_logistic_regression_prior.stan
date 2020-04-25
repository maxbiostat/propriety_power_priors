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
transformed parameters{
  real logL = bernoulli_logit_lpmf(y0 | alpha + X0*beta);
  real logL_sq = square(logL);
}
model {
  /*power prior*/
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta | 0, 1);
  target += a_0 * logL;
}
