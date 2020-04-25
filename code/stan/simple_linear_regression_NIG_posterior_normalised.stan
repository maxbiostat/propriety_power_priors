data{
  int<lower=0> P;
  int<lower=0> N0;
  vector[N0] y0;
  matrix[N0, P] X0;
  vector[P] mu_beta;
  matrix[P, P] lambda_0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real<lower=0> nu;
  real<lower=0> eta;
  int<lower=0> N;
  real y[N];
  matrix[N, P] X;
}
transformed data{
  matrix[P, P] invlambda_0 = inverse_spd(lambda_0);
}
parameters{
  vector[P] beta;
  real<lower=0> sigma_sq;
  real<lower=0, upper=1> a_0;
}
model{
  matrix[N0, P] Xstar = sqrt(a_0) * X0;
  vector[N0] ystar = sqrt(a_0) * y0;
  matrix[P, P] lambda_n = (Xstar'*Xstar) + invlambda_0;   
  vector[P] mu_n = inverse_spd(lambda_n) * (invlambda_0*mu_beta + Xstar'*ystar);
  real alpha_n =  alpha0 + (N0*a_0)/2.0;
  real beta_n = beta0 + 0.5 * ( ystar'*ystar + mu_beta'*invlambda_0*mu_beta - mu_n'*lambda_n*mu_n );
  real lconst = -(a_0 * N0/2) * log(2*pi()) + .5 * (log(determinant(invlambda_0)) - log(determinant(lambda_n))) + (alpha0 * log(beta0) - alpha_n * log(beta_n) ) + (lgamma(alpha_n)-lgamma(alpha0));
  /*power prior*/
  target += a_0 * normal_lpdf(y0 | X0*beta, sqrt(sigma_sq) );
  target += multi_normal_lpdf(beta| mu_beta, sigma_sq * lambda_0);
  target += inv_gamma_lpdf(sigma_sq | alpha0, beta0);
  target += -lconst;
  target += beta_lpdf(a_0 | eta, nu);
  /* Likelihood */
  target += normal_lpdf(y | X*beta, sqrt(sigma_sq) );
}
