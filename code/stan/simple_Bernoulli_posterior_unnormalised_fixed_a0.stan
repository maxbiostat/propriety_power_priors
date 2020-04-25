data{
  int<lower=1> N0;
  int<lower=0> y0;
  real<lower=0> c;
  real<lower=0> d;
  int<lower=1> N;
  int<lower=0> y;
  real<lower=0> a_0;
}
parameters{
  real<lower=0,upper=1> theta;
}
transformed parameters{
  real logL = y0 * log(theta) + (N0 - y0)  * log1m(theta);
  real logL_sq = square(logL);
}
model{
  /*Power prior*/
  target += a_0 * logL;   
  target += beta_lpdf(theta | c, d);
  /* Likelihood */
  target += y * log(theta) + (N - y)  * log1m(theta);
}
