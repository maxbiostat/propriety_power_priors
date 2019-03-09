data{
  int<lower=1> N0;
  int<lower=0> y0;
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0> delta;
  real<lower=0> nu;
  int<lower=1> N;
  int<lower=0> y;
}
parameters{
  real<lower=0,upper=1> theta;
  real<lower=0,upper=1> a_0;
}
model{
  /*Power prior*/
  target += a_0 * binomial_lpmf(y0 | N0, theta);
  target += beta_lpdf(theta | c, d);
  target += beta_lpdf(a_0 | delta, nu);
  /* Likelihood */
  target += binomial_lpmf(y | N, theta);
}
