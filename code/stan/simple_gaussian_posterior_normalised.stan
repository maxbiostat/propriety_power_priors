data{
  int<lower=1> N0;
  real y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real mu0;
  real<lower=0> kappa0;
  real<lower=0> delta;
  real<lower=0> nu;
  int<lower=0> N;
  real y[N];
}
transformed data{
  real ybar;
  real s = 0;
  ybar = mean(y0);
  for(i in 1:N0) s += (1.0/N0) * (y0[i] - ybar)^2;
}
parameters{
  real mu;
  real<lower=0> sigma_sq;
  real<lower=0,upper=1> a_0;
}
model{
  real nstar = a_0 * N0;
  real alpha_n = alpha0 + nstar/2;
  real kappa_n = kappa0 + nstar;
  real beta_n = beta0 + .5 * (nstar * s +  (kappa0 * nstar * (ybar - mu0)^2 )/kappa_n );
  real term1 = lgamma(alpha_n)-lgamma(alpha0);
  real term2 = alpha0 * log(beta0) - alpha_n * log(beta_n) ;
  real term3 = 0.5 *( log(kappa0) - log(kappa_n) )-nstar/2 * log(2*pi()); 
  /* Power prior */
  target += a_0 * normal_lpdf(y0 | mu, sqrt(sigma_sq));
  target += normal_lpdf(mu | mu0, sqrt( 1/kappa0 * sigma_sq ) );
  target += inv_gamma_lpdf(sigma_sq | alpha0, beta0);
  target += beta_lpdf(a_0 | delta, nu);
  // normalising constant stuff
  target += -( term1 + term2 + term3);
  /* Likelihood */
  target += normal_lpdf(y | mu, sqrt(sigma_sq));
}
