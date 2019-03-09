data{
  int<lower=1> N0;
  real y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real mu0;
  real<lower=0> kappa0;
  real<lower=0> a_0;
}
parameters{
  real mu;
  real<lower=0> sigma_sq;
}
model{
  target += a_0 * normal_lpdf(y0 | mu, sqrt(sigma_sq));
  target += normal_lpdf(mu | mu0, sqrt( 1/kappa0 * sigma_sq ) );
  target += inv_gamma_lpdf(sigma_sq | alpha0, beta0);
}
