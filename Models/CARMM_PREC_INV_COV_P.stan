data {
  // Number of areas
  int<lower=0> n;
  
  // Number of memberships
  int<lower = 0> m;

  // Spatial matrix
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  
  // Covariates
  int<lower=0> k;   // number of predictors
  matrix[m, k] X_cov;   // predictor matrix
  
  // Outcomes
  int<lower = 0> y[m];
  
  // Offsets
  vector[m] log_offset;
  
  // Inverse Membership
  matrix[n, m] I_H;
}
transformed data {
  // For implied adjacency by MM
  matrix[m, m] W_t;
  matrix[m, m] D_t;
  
  // Adjacency matrix
  matrix<lower = 0>[n, n] D;
  {
    vector[n] W_rowsums;
    for (i in 1:n) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
  
  // Implied adjacency by MM
  D_t = I_H' * D * I_H;
  W_t = I_H' * W * I_H;
}
parameters {
  real<lower = 0, upper = 1> alpha;
  real<lower = 0> tau;
  vector[m] phi_unsc;
  
  // Intercept and covariates
  real gamma;
  vector[k] beta; // coefficients on Q_ast
}
transformed parameters{
  vector[m] r_1;
  real<lower = 0> invtausq;
  vector[m] phi_mm;
  
  invtausq = inv_sqrt(tau);
  phi_mm = invtausq*phi_unsc;
  
  // Relative risk
  r_1 = gamma + X_cov * beta + phi_mm;
}
model {
  // CAR
  phi_unsc ~ multi_normal_prec(rep_vector(0, m), (D_t - alpha * W_t));
  tau ~ gamma(2, .2);
  
  // Linear term
  beta ~ normal(0, .7);
  gamma ~ normal(0, .7);
  
  // Poisson
  y ~ poisson(exp(log_offset + r_1));
}
generated quantities {
  int<lower = 0> yrep[m];
  vector[m] log_lik;
  vector[m] log_lik_rep;
  real sum_ll;
  real sum_ll_rep;
  real ppp;
  vector[n] l_RR;
  vector[n] phi;
  
  // Retrieve areal relative risk
  l_RR = I_H * r_1;
  phi = I_H * phi_mm;
  
  // Simulate from posterior
  for (i in 1:m) {
  //   // likelihood of the current parameter values (given the original data)
    log_lik[i] = poisson_lpmf(y[i] | exp(r_1[i] + log_offset[i]));
  //   // generate new data based on current parameter values
    yrep[i] = poisson_rng(exp(r_1[i] + log_offset[i]));
  //   // compute likelihood of the current parameter values (given the new data)
    log_lik_rep[i] = poisson_lpmf(yrep[i] | exp(r_1[i] + log_offset[i]));
  }
  //   // sum up the likelihoods for all observations
    sum_ll = sum(log_lik) ;
    sum_ll_rep = sum(log_lik_rep);
  //   // check which is higher
    ppp = sum_ll > sum_ll_rep ? 1 : 0;
}
