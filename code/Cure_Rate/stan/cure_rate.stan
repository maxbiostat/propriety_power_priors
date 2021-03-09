functions{
  real Surv(real t, real alpha, real lambda){
    return(exp(-exp(lambda) * t^alpha));
  }
  real cure_rate_likelihood_lpdf(real x, real delta, real alpha, real lambda, real theta){
    // Latent N's summed out
    real ans;
    real S = Surv(x, alpha, lambda);
    real logF = weibull_lpdf(x | alpha, exp(-lambda/alpha)) ;
    ans = delta*(log(theta) + logF) - theta * (1 - S) ;
    return(ans); 
  }
}
data{
  int<lower=0> n;
  int<lower=0> p;
  real<lower=0> y[n];
  matrix[n, p] X;
  int<lower=0, upper=1> delta[n];
  /*Hyper parameters*/
  real mu_0;
  real sigma_0;
  real<lower=0> delta_0;
  real<lower=0> tau_0;
  real<lower=0> sigma_beta;
}
parameters{
  real lambda;
  real<lower=0> alpha;
  vector[p] beta;
}
transformed parameters{
  real logL = 0;
  real logL_sq = 0;
  {
    vector[n] theta = exp(X*beta);
    for(i in 1:n){
      logL += cure_rate_likelihood_lpdf(y[i] | delta[i], alpha, lambda, theta[i]);
    }
  }
  logL_sq = square(logL);
}
model{
  target += logL;
  /*Priors*/
  target += normal_lpdf(lambda | mu_0, sigma_0);
  target += gamma_lpdf(alpha | delta_0, tau_0);
  target += normal_lpdf(beta | 0, sigma_beta);
}
