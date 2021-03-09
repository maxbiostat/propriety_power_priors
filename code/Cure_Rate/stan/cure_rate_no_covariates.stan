functions{
  real Surv(real t, real alpha, real sigma){
    return(exp(weibull_lccdf(t  | alpha, sigma)));
  }
  real cure_rate_likelihood_lpdf(real x, real delta, real alpha, real lambda, real theta){
    // Latent N's summed out
    real ans;
    real sigma = exp(-lambda/alpha);
    real logF = weibull_lpdf(x | alpha, sigma) ;
    real S = Surv(x, alpha, sigma);
   ans = delta*(log(theta) + logF) - theta * (1 - S) ;
   return(ans); 
  }
}
data{
  int<lower=0> n;
  real<lower=0> y[n];
  int<lower=0, upper=1> delta[n];
  /*Hyper parameters*/
  real mu_0;
  real sigma_0;
  real<lower=0> delta_0;
  real<lower=0> tau_0;
}
parameters{
  real lambda;
  real<lower=0> alpha;
  real<lower=0> theta;
}
model{
  for(i in 1:n){
    target += cure_rate_likelihood_lpdf(y[i] | delta[i], alpha, lambda, theta);
  }
  /*Priors*/
  target += normal_lpdf(lambda | mu_0, sigma_0);
  target += gamma_lpdf(alpha | delta_0, tau_0);
  /* Temporary */
  target += normal_lpdf(theta | 0, 100);
}
