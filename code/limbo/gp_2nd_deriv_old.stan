functions {
  real K22(real t, real tp, real alpha, real rho){
    real cov; // computation in log-space for numerical stability
    real sq_delta = square(t -tp);
    real sq_rho = square(rho);
    /* potentially numerically unstable...*/
    cov = square(alpha) * exp(-sq_delta/(2*sq_rho)) * 1/square(sq_rho) * ( square(sq_delta)/square(sq_rho) - (6*sq_delta/sq_rho) + 3 );
    return(cov);
  }
  matrix K22_cov_mat(real[] s,  real alpha, real rho, int K){
    matrix[K, K] M;
    real sq_diag = 3*square(alpha/square(rho));
    for (i in 1:K) {
      M[i, i] = sq_diag;
      for (j in 1:K) {
        M[i, j] = K22(s[i], s[j], alpha, rho);
        M[j, i] = M[i, j];
      }
    }
     return(cholesky_decompose(M));
  }
  vector gp_pred_rng(real[] x2,
                     vector y1, real[] x1,
                     real alpha, real rho, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K = cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
    vector gp_mu(real[] x2,
                     vector y1, real[] x1,
                     real alpha, real rho, real sigma) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
      matrix[N1, N1] K = cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      return(f2_mu);
  }
}
data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
  int<lower=1> K;
  real s[K];
}

transformed data {
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector<lower=0>[K] zpp; // Second derivative process
}
transformed parameters{
}

model {    
  matrix[N, N] cov = cov_exp_quad(x, alpha, rho) + diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
  matrix[K, K] cM = K22_cov_mat(s, alpha, rho, K);
  y ~ multi_normal_cholesky(rep_vector(0, N), L_cov);
  /* Priors*/
  rho ~ normal(0, 1);
  alpha ~ normal(0, 10);
  sigma ~ normal(0, 1);
  zpp ~ multi_normal_cholesky(rep_vector(0, K), cM);
}
generated quantities{
  vector[K] f_predict = gp_pred_rng(s, y, x, alpha, rho, sigma, 1e-10);
  vector[K] f_mu = gp_mu(s, y, x, alpha, rho, sigma);
  vector[K] y_predict;
  for (n in 1:K)
    y_predict[n] = normal_rng(f_predict[n], sigma);
}
