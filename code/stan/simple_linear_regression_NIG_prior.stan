data{
  int<lower=0> N0;
  int<lower=0> P;
  real y0[N0];
  matrix[N0, P] X0;
  vector[P] mu_beta;
  matrix[P, P] lambda_0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real<lower=0> a_0;
}
parameters{
  vector[P] beta;
  real<lower=0> sigma_sq;
}
transformed parameters{
  real logL = normal_lpdf(y0 | X0*beta, sqrt(sigma_sq) );
  real logL_sq = square(logL);
}
model{
  /*power prior*/
  target += multi_normal_lpdf(beta| mu_beta, sigma_sq * lambda_0);
  target += inv_gamma_lpdf(sigma_sq | alpha0, beta0);
  target += a_0 * logL;
}
