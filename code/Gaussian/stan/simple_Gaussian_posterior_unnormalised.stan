data{
  int<lower=1> N0;
  real y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real mu0;
  real<lower=0> kappa0;
  real<lower=0> nu;
  real<lower=0> eta;
  int<lower=0> N;
  real y[N];
}
parameters{
  real mu;
  real<lower=0> tau;
  real<lower=0, upper=1> a_0;
}
transformed parameters{
  real logL = normal_lpdf(y0 | mu, sqrt(1/tau));
  real logL_sq = square(logL);
}
model{
  /* Power prior */
  target += a_0 * logL;
  target += normal_lpdf(mu | mu0, sqrt( 1/(kappa0 * tau)  ) );
  target += gamma_lpdf(tau | alpha0, beta0);
  target += beta_lpdf(a_0 | eta, nu);
  /* Likelihood */
  target += normal_lpdf(y | mu, sqrt(1/tau));
}
generated quantities{
  real<lower=0> sigma = 1/tau;
}
