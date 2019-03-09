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
  real<lower=0> a_0;
}
parameters{
  real<lower=0> mu;
  real<lower=0> lambda;
}
model{
  /*Power prior*/
  for( i in 1:N0) target += a_0 * inverse_gaussian_lpdf(y0[i] | mu, lambda);
  target += gamma_lpdf(mu | alpha_m, beta_m);
  target += gamma_lpdf(lambda | alpha_l, beta_l);
}
