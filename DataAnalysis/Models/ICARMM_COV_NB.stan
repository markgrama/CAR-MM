functions {
  real sparse_iar_lpdf(vector phi, 
  real tau,
  // real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      // for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      // return 0.5 * (n * log(tau)
      //               + sum(ldet_terms)
      //               - tau * (phit_D * phi - alpha * (phit_W * phi)));
      return 0.5 * ((n-1) * log(tau) - tau * (phit_D * phi - (phit_W * phi)));
  }
  real corr_val(vector x, vector y, int n){
    int n_vec;
    real num;
    real den_1;
    real den_2;
    n_vec = n;
    num = n_vec*sum(x .* y) - sum(x)*sum(y);
    den_1 = sqrt(n_vec*sum(x .* x) - pow(sum(x), 2));
    den_2 = sqrt(n_vec*sum(y .* y) - pow(sum(y), 2));
    return num/(den_1 * den_2);
    }
}
data {
  // Number of areas
  int<lower=0> n;

  // Multiple membership
  int<lower=0> m;
  matrix[m, n] H;
  
  // Spatial matrix
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  int W_n; // number of adjacent region pairs
  
  // Covariates
  int<lower=0> k;   // number of predictors
  matrix[n, k] X_cov;   // predictor matrix
  
  // Outcomes
  int<lower = 0> y[m];
  
  // Offsets
  vector[m] log_offset;
}
transformed data {
   // Adjacency
  // int W_n; // number of adjacent region pairs
  // Number of directed neighbours per area
  vector[n] n_i;
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:n) D_sparse[i] = sum(W[i]);
  {
    vector[n] invsqrtD;  
    for (i in 1:n) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}
parameters {
  real<lower = 0> tau;
  vector[n] phi1;
  
  // Intercept and covariates
  real gamma;
  vector[k] beta;
  // Negative Binomial
  real<lower = 0> psi;
}
transformed parameters{
  vector[m] r_1;
  real<lower = 0> invtausq;
  vector[n] phi;

  invtausq = inv_sqrt(tau);
  phi = invtausq*phi1;
  
  // Relative risk
  r_1 = H*(gamma + X_cov * beta + phi);
}
model {
  // y ~ poisson_log(log_E + beta0 + beta1 * x + phi * sigma_phi + theta * sigma_theta);

  // NOTE:  no prior on phi_raw, it is used to construct phi
  // the following computes the prior on phi on the unit scale with sd = 1
  target += -0.5 * dot_self(phi1[W_sparse[,1]] - phi1[W_sparse[,2]]);
  // soft sum-to-zero constraint on phi)
  sum(phi1) ~ normal(0, 0.001 * n);  // equivalent to mean(phi) ~ normal(0,0.001)
  
  psi ~ gamma(2, .2);
  tau ~ gamma(2, .2);
  // psi ~ normal(0, 10);
  // tau ~ normal(0, 10);
  beta ~ normal(0, .7);
  gamma ~ normal(0, .7);
  
  y ~ neg_binomial_2(exp(log_offset + r_1), psi);
}
generated quantities {
  int<lower = 0> yrep[m];
  vector[m] log_lik;
  vector[m] log_lik_rep;
  real sum_ll;
  real sum_ll_rep;
  real ppp;
  vector[n] l_RR;
  
  // Areal relative risk
  l_RR = gamma + X_cov * beta + phi;
  
  // Simulate from posterior
  for (i in 1:m) {
    // likelihood of the current parameter values (given the original data)
    log_lik[i] = neg_binomial_2_lpmf(
      y[i] | exp(r_1[i] + log_offset[i]), psi
    );
    // generate new data based on current parameter values
    yrep[i] = neg_binomial_2_rng(
      exp(r_1[i] + log_offset[i]), psi
    );
    // compute likelihood of the current parameter values (given the new data)
    log_lik_rep[i] = neg_binomial_2_lpmf(
      yrep[i] | exp(r_1[i] + log_offset[i]), psi
    );
  }
    // sum up the likelihoods for all observations
    sum_ll = sum(log_lik) ;
    sum_ll_rep = sum(log_lik_rep);
    // check which is higher
    ppp = sum_ll > sum_ll_rep ? 1 : 0;
}
