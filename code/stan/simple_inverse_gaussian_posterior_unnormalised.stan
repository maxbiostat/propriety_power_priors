functions{
  real inverse_gaussian_lpdf(real x, real mu, real lambda){
    real ldens = 0.5 * ( log(lambda) - log(2*pi()) - 3*log(x) ) - (lambda *(x-mu)^2)/(2* mu^2 * x); 
    return(ldens);
  }
}
data{
  int<lower=0> N0;
  real<lower=0> y0[N0];
  real<lower=0> alpha_l;
  real<lower=0> beta_l;
  real<lower=0> alpha_m;
  real<lower=0> beta_m;
  int<lower=0> N;
  real<lower=0> y[N];
  real<lower=0> delta;
  real<lower=0> nu;
}
parameters{
  real<lower=0> mu;
  real<lower=0> lambda;
  real<lower=0,upper=1> a_0;
}
transformed parameters{
  real logL = inverse_gaussian_lpdf(y0[i] | mu, lambda);
  real logL_sq = square(logL);
}
model{
  /*Power prior*/
  for(i in 1:N0) target += a_0 * logL;
  target += gamma_lpdf(mu | alpha_m, beta_m);
  target += gamma_lpdf(lambda | alpha_l, beta_l);
  target += beta_lpdf(a_0 | delta, nu);
  /* Likelihood */
  for(i in 1:N) target += inverse_gaussian_lpdf(y[i] | mu, lambda);
}
