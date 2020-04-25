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
transformed parameters{
  real logL = 0;
  real logL_sq = 0;
  for(i in 1:N0) logL += inverse_gaussian_lpdf(y0[i] | mu, lambda);
  logL_sq = square(logL);
}
model{
  /*power prior*/
  target += gamma_lpdf(mu | alpha_m, beta_m);
  target += gamma_lpdf(lambda | alpha_l, beta_l);
  target += a_0*logL;
}
