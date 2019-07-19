functions {
  /*Functions to create all three covariance (sub) matrices (fourth is just cov_quad_exp) */
  real K22(real s, real sp, real alpha, real rho){
    real sq_delta = square(s-sp);
    real sq_rho = square(rho);
    real cov = square(alpha) * exp(-sq_delta/(2*sq_rho)) * 1/square(sq_rho) * ( square(sq_delta)/square(sq_rho) - (6*sq_delta/sq_rho) + 3 );
    return(cov);
  }
  matrix K22_cov_mat(real[] s,  real alpha, real rho){
    int K = size(s);
    matrix[K, K] X;
    real diag_entry = 3*square(alpha/square(rho));
    for (i in 1:K) {
      X[i, i] = diag_entry;
      for (j in 1:K) {
        X[i, j] = K22(s[i], s[j], alpha, rho);
        X[j, i] = X[i, j];
      }
    }
    return(X);
  }
  real K20(real t, real s, real alpha, real rho){
    real cov; // computation in log-space for numerical stability
    real sq_delta = square(t-s);
    real sq_rho = square(rho);
    cov = square(alpha) *  exp(-sq_delta/(2*sq_rho)) * ( sq_delta / square(sq_rho) - 1/sq_rho);
    return(cov);
  }
  matrix K20_cov_mat(real[] x, real[] s, real alpha, real rho){
    int N = size(x);
    int K = size(s);
    matrix[N, K] Y;
    for (i in 1:N) {
      for (j in 1:K) {
        Y[i, j] = K20(x[i], s[j], alpha, rho);
      }
    }
    return(Y);
  }
  matrix full_cov_mat(real [] x, real[] s, real alpha, real rho, real epsilon){
    int N = size(x);
    int K = size(s);
    int J = N + K;
    matrix[N, N] A;
    matrix[N, K] B;
    matrix[K, K] C;
    matrix[J, J] M;
    // | A, B|
    // | B',C|
    A = cov_exp_quad(x, alpha, rho);
    B = K20_cov_mat(x, s, alpha, rho);
    C = K22_cov_mat(s, alpha, rho);
    /*Filling in the matrix*/
    M[1:N, 1:N] = A;
    M[1:N, (N+1):J] = B;
    M[(N + 1):J, (N+1):J] = C;
    M[(N + 1):J, 1:N] = B';
    return(M + diag_matrix(rep_vector(epsilon,  J)) );
  }
  /* Prediction stuff*/
  vector gp_pred_rng(real[] x2,
  vector y1,
  real[] x1,
  real alpha,
  real rho,
  real sigma,
  real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] L_K;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      K = cov_exp_quad(x1, alpha, rho);
      for (n in 1:N1)
      K[n, n] = K[n,n] + square(sigma);
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      f2_mu = (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      cov_f2 = cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}
data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
  int<lower=1> K;
  real s[K];
  int m[K];
  real<lower=0> epsilon;
  int<lower=1> N_pred;
  real x_pred[N_pred];
}
transformed data {
  int N_tot = N + K;
  vector[N_tot] mu = rep_vector(0, N_tot);
  real delta_x = max(x)-min(x);
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] f_obs;
  vector<lower=0>[K] fpp;
}
transformed parameters{
  vector[N_tot] f;
  vector[N_tot] f_tilde;
  f_tilde[1:N] = f_obs;
  f_tilde[(N+1):N_tot] = fpp;
  f =  cholesky_decompose(full_cov_mat(x, s, alpha, rho, epsilon)) * f_tilde;
}

model {
  y ~ normal(f[1:N], sigma);
  /* Priors*/
  rho ~ normal(0, delta_x/3.0);
  alpha ~ normal(0, 1);
  sigma ~ cauchy(0, 10);
  f_tilde ~ normal(0, 1);
  f[1] ~ normal(0, 1E-16);
}
generated quantities{
  vector[N_pred] f_pred = gp_pred_rng(x_pred, y, x, alpha, rho, sigma, epsilon);
}
