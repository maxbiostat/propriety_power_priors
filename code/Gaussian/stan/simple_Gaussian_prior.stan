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
  real<lower=0> tau;
}
transformed parameters{
  real logL =  normal_lpdf(y0 | mu, sqrt(1/tau) );
  real logL_sq = square(logL);
}
model{
  /*power prior*/
  target += normal_lpdf(mu | mu0, sqrt( 1/(kappa0 * tau)  ) );
  target += gamma_lpdf(tau | alpha0, beta0);
  target += a_0 * logL;
}
generated quantities{
  real<lower=0> sigma = sqrt(1/tau);
}
