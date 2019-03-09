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
model{
  target += a_0 * poisson_lpmf(y0 | lambda) ;
  target += gamma_lpdf(lambda | alpha0, beta0);
}
