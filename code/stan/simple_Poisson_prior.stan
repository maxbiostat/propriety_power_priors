data{
  int<lower=0> N0;
  int y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real<lower=0> a_0;
}
parameters{
  real<lower=0> lambda;
}
transformed parameters{
  real logL = poisson_lpmf(y0 | lambda);
  real logL_sq = square(logL);
}
model{
  /*power prior*/
  target += gamma_lpdf(lambda | alpha0, beta0);
  target += a_0 * logL;
}
