#### SETUP ####

# DATASET INDEX (WHICH OF THE K DATASETS TO SAMPLE FROM)
index <- 1
# SIMULATION IDENTIFIER (GENERAL ROUND OF SIMULATIONS)
sim_id <- "check"
# HMC PARAMETERS
niter <- 1E3
nchains <- 3
# POST or INVERSE DATA GENERATION PARAMETERISATION
param <- "post"
# DIMENSION
dim_par <- 100
# MISSING TRIGGER
missing <- T
# THINNING PARAMETER
par_thinning <- 5

# TODO: SET WORKING DIRECTORY
setwd(path)
source("Functions.R")

# LOAD DATA ACCORDING TO PARAMETERISATION
if(param == "post"){
  load(file.path(path, "Data/Data_INV_POST_SBC_POST.Rdata"))
  sim_id <- paste(sim_id, "_", dim_par, "_post", sep = "")
} else{
  if(param == "inv"){
    load(file.path(path, "HPC/INV_POST_SBC/Data/Data_INV_POST_SBC_INV.Rdata"))
    sim_id <- paste(sim_id, "_", dim_par, "_inv", sep = "")
  }
}

# Libraries
library(rstan)
## CHOOSE WHETHER TO SET auto_write = T/F BASED ON YOUR SETUP
library(dplyr)
## PARALLEL OPTIONS IF NEEDED
# options(mc.cores = ncores)
options(mc.cores = parallel::detectCores())
print(paste("ID:", sim_id))
#

#### SAMPLING ####

# Fetch dataset
d_grid <- get(paste("d_grid_", dim_par, sep = ""))
## INVERSE
if(dim_par != 130){
  d_inv <- list(
    # Number of areas
    n = nrow(d_grid$W),
    # Multiple membership
    m = nrow(d_grid$H),
    I_H = d_grid$H_inv,
    # Covariates
    k = ncol(d_grid$X_mm),
    X_cov = d_grid$X_mm,
    # Adjacency
    W_n = sum(d_grid$W) / 2,
    W = d_grid$W,
    # Outcome
    y = d_grid$y[index,],
    log_offset = log(d_grid$off_grid)
  )
  ## MODEL
  m_inv <- stan_model(
    file = "Models/CARMM_PREC_INV_COV_P.stan"
  )
  ## SAMPLING
  time_inv <- Sys.time()
  fit_inv <- sampling(
    m_inv, data = d_inv,
    iter = niter, chains = nchains,
    control = list(adapt_delta = .99, max_treedepth = 15)
  )
  time_inv <- Sys.time() - time_inv
}
## POST
d_post <- list(
  # Number of areas
  n = nrow(d_grid$W),
  # Multiple membership
  m = nrow(d_grid$H),
  H = d_grid$H,
  # Covariates
  k = ncol(d_grid$X_cov),
  X_cov = d_grid$X_cov,
  # Adjacency
  W_n = sum(d_grid$W) / 2,
  W = d_grid$W,
  # Outcome
  y = d_grid$y[index,],
  log_offset = log(d_grid$off_grid)
)

# Compile models
m_post <- stan_model(file = "Models/CARMM_COV_P.stan")
## POST
time_post <- Sys.time()
fit_post <- sampling(
  m_post, data = d_post,
  iter = niter, chains = nchains,
  control = list(adapt_delta = .99, max_treedepth = 15)
)
time_post <- Sys.time() - time_post
#

#### GET RESULTS ####

# Global parameters
## General dimentions
B <- niter*nchains/2 ; n <- nrow(d_grid$W) ; m <- nrow(d_grid$H)
## Names
### CAR
phi_names <- sapply(1:n, function(x) paste("phi[", x, "]", sep = ""))
### Areal RR
l_RR_names <- sapply(1:n, function(x) paste("l_RR[", x, "]", sep = ""))
RR_names <- sapply(1:n, function(x) paste("RR[", x, "]", sep = ""))
### Membership RR
rr_mm_names <- sapply(1:m, function(x) paste("r_1[", x, "]", sep = ""))
### Mean function parameters
pars_mean <- c("alpha", "tau", "gamma", "beta[1]", "beta[2]")
### Final vector
all_pars <- c(pars_mean, phi_names, l_RR_names, rr_mm_names)

if(dim_par != 130){
  post_samples_inv <- f_ar_mat(
    fit = fit_inv, inv = F, niter = niter, nchains = nchains, 
    pars_mean = pars_mean, phi_names = phi_names, l_RR_names = l_RR_names,
    rr_mm_names = rr_mm_names
  )
  array_inv <- post_samples_inv$array ; mat_inv <- post_samples_inv$mat
  
  # SBC
  ## INVERSE
  results_inv <- f_sbc(
    fit = fit_inv, inv = F, niter = niter, nchains = nchains, 
    pars_mean = pars_mean, phi_names = phi_names, l_RR_names = l_RR_names,
    RR_names = RR_names, rr_mm_names = rr_mm_names, list_pars = d_grid,
    index = index, array_fit = array_inv, mat_fit = mat_inv,
    thinning = par_thinning
  )
} else{
  results_inv <- NULL ; time_inv <- NULL
}

## POST
post_samples_post <- f_ar_mat(
  fit = fit_post, inv = F, niter = niter, nchains = nchains, 
  pars_mean = pars_mean, phi_names = phi_names, l_RR_names = l_RR_names,
  rr_mm_names = rr_mm_names
)
array_post <- post_samples_post$array ; mat_post <- post_samples_post$mat

## POST
results_post <- f_sbc(
  fit = fit_post, inv = F, niter = niter, nchains = nchains, 
  pars_mean = pars_mean, phi_names = phi_names, l_RR_names = l_RR_names,
  RR_names = RR_names, rr_mm_names = rr_mm_names, list_pars = d_grid,
  index = index, array_fit = array_post, mat_fit = mat_post,
  thinning = par_thinning
)
#

#### FINAL SAVE ####

# Results
res_folder <- file.path(
  "Results/Posterior_samples"
)

# Iterations for SBC plots
l_iter <- list(niter = niter, nchains = nchains, par_thinning = par_thinning)

save(
  results_inv, results_post, l_iter,
  time_inv, time_post,
  file = file.path(
    res_folder, paste("sbc_", index, "_", sim_id, ".RData", sep = "")
  )
)
#




