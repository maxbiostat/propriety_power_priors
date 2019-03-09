data{
  int<lower=1> N0;
  int<lower=0> y0;
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0> a_0;
}
parameters{
  real<lower=0,upper=1> theta;
}
model{
  target += a_0 * binomial_lpmf(y0 | N0, theta);
  target += beta_lpdf(theta | c, d);
}
