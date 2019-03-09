data{
  int<lower=0> P;
  int<lower=0> N0;
  real y0[N0];
  matrix[N0, P] X0;
  vector[P] mu_beta;
  matrix[P, P] lambda_0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real<lower=0> nu;
  real<lower=0> delta;
  int<lower=0> N;
  real y[N];
  matrix[N, P] X;
}
parameters{
  vector[P] beta;
  real<lower=0> sigma_sq;
  real<lower=0,upper=1> a_0;
}
model{
  /*power prior*/
  target += a_0 * normal_lpdf(y0 | X0*beta, sqrt(sigma_sq) );
  target += multi_normal_lpdf(beta| mu_beta, sigma_sq * lambda_0);
  target += inv_gamma_lpdf(sigma_sq | alpha0, beta0);
  target += beta_lpdf(a_0 | nu, delta);
  /* Likelihood */
  target += normal_lpdf(y | X*beta, sqrt(sigma_sq) );
}
