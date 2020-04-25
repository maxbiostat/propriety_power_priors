data{
  int<lower=0> N0;
  int y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  int<lower=0> N;
  int y[N];
  real<lower=0> nu;
  real<lower=0> eta;
}
transformed data{
  real s = sum(y0);
  real logPprime = 0.0;
  for(i in 1:N0) logPprime += lgamma(y0[i] + 1);
}
parameters{
  real<lower=0> lambda;
  real<lower=0, upper=1> a_0;
}
model{
  /* Power prior */
  target += -( -a_0 * logPprime + lgamma(a_0 * s + alpha0) - (a_0 * s + alpha0)* log(a_0 * N0 + beta0) );
  target += a_0 * poisson_lpmf(y0 | lambda) ;
  target += gamma_lpdf(lambda | alpha0, beta0);
  target += beta_lpdf(a_0 | eta, nu);
  /* Likelihood */
  target += poisson_lpmf(y | lambda);
}
