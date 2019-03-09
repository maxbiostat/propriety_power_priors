data{
  int<lower=0> N0;
  real<lower=0> y0[N0];
  real<lower=0> eta_a;
  real<lower=0> nu_a;
  real<lower=0> eta_b;
  real<lower=0> nu_b;
  real<lower=0> a_0;
}
parameters{
  real<lower=0> alpha;
  real<lower=0> beta;
}
model{
  target += a_0 * gamma_lpdf(y0 | alpha, beta);
  target += gamma_lpdf(alpha | eta_a, nu_a);
  target += gamma_lpdf(beta | eta_b, nu_b);
}
