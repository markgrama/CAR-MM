data {
  // Number of areas
  int<lower=0> n;

  // Multiple membershiop
  int<lower=0> m;
  matrix[m, n] H;

  // Covariates
  int<lower=0> k;   // number of predictors
  matrix[n, k] X_cov;   // predictor matrix
  
  // Outcomes
  int<lower = 0> y[m];
  
  // Offsets
  vector[m] log_offset;
}
parameters {
  // Negative Binomial overdispersion
  real<lower = 0> psi;
  
  // Intercept and covariates
  real gamma;
  vector[k] beta;
}
transformed parameters{
  vector[m] r_mm;
  
  // Relative risk
  r_mm = H*(gamma + X_cov * beta);
}
model {
  beta ~ normal(0, .7);
  gamma ~ normal(0, .7);
  
  // Negative Binomial
  psi ~ gamma(2, .2);
  y ~ neg_binomial_2(exp(log_offset + r_mm), psi);
}
generated quantities {
  int<lower = 0> yrep[m];
  vector[m] log_lik;
  vector[m] log_lik_rep1;
  real sum_ll;
  real sum_ll_rep;
  real ppp;
  vector[n] l_RR;
  
  // Areal relative risk
  l_RR = gamma + X_cov * beta;
  
  // Simulate from posterior
  for (i in 1:m) {
    // likelihood of the current parameter values (given the original data)
    log_lik[i] = neg_binomial_2_lpmf(y[i] | exp(r_mm[i] + log_offset[i]), psi);
    // generate new data based on current parameter values
    yrep[i] = neg_binomial_2_rng(exp(r_mm[i] + log_offset[i]), psi);
    // compute likelihood of the current parameter values (given the new data)
    log_lik_rep1[i] = neg_binomial_2_lpmf(yrep[i] | exp(r_mm[i] + log_offset[i]),
                                                    psi);
  }
    // sum up the likelihoods for all observations
    sum_ll = sum(log_lik) ;
    sum_ll_rep = sum(log_lik_rep1);
    // check which is higher
    ppp = sum_ll > sum_ll_rep ? 1 : 0;
}
