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
  int<lower=0> P;
  int<lower=0> N0;
  matrix[N0, P] X0;
  real<lower=0> y0[N0];
  int<lower=0, upper=1> Delta0[N0];
  int<lower=0> N;
  matrix[N, P] X;
  real<lower=0> y[N];
  int<lower=0, upper=1> Delta[N];
  /*Hyper parameters*/
  real mu_0;
  real sigma_0;
  real<lower=0> delta_0;
  real<lower=0> tau_0;
  real<lower=0> sigma_beta;
  real<lower=0> eta;
  real<lower=0> nu;
}
parameters{
  real lambda;
  real<lower=0> alpha;
  vector[P] beta;
  real<lower=0, upper=1> a_0;
}
transformed parameters{
  real logL0 = 0;
  real logL = 0;
  {
    vector[N0] theta0 = exp(X0*beta);
    vector[N] theta = exp(X*beta);
      for(j in 1:N0) logL0 += cure_rate_likelihood_lpdf(y0[j] | Delta0[j], alpha, lambda, theta0[j]);
      for(i in 1:N) logL += cure_rate_likelihood_lpdf(y[i] | Delta[i], alpha, lambda, theta[i]);
  }
}
model{
  target += logL;
  /*Priors*/
  target += a_0 * logL0;
  target += normal_lpdf(lambda | mu_0, sigma_0);
  target += gamma_lpdf(alpha | delta_0, tau_0);
  target += normal_lpdf(beta | 0, sigma_beta);
  target += beta_lpdf(a_0 | eta, nu);
}
