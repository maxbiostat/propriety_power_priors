data {
  int<lower=0> D;
  int<lower=0> N0max;
  int<lower=0> N0[D]; // where each array ends
  int<lower=0> P;
  real X0[N0max, P, D];
  int<lower=0, upper=1> y0[N0max, D];
  int<lower=0> N;
  matrix[N, P] X;
  int<lower=0, upper=1> y[N];
  real<lower=0> eta[D];
  real<lower=0> nu[D];
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0,upper=1> a_0[D];
}
model {
  /* Power prior */
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta | 0, 1);
  for(d in 1:D){
    target += a_0[d] * bernoulli_logit_lpmf( y0[1:N0[d], d] | alpha + to_matrix(X0[1:N0[d], , d])*beta);
    target += beta_lpdf(a_0[d] | eta[d], nu[d]);
  }  
  /* Likelihood */
  target += bernoulli_logit_lpmf(y | alpha + X*beta);
}
