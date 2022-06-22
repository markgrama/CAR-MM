#### GENERAL SETUP ####

# TODO: SET WORKING DIRECTORY
setwd(file.path(path))

# Load package
## Data handling
library(dplyr)
## Sampling
library(MASS)
## Matrix
library(Matrix)

# Load functions
source("Functions.R")
#

#### COMMON ####

# General parameters
K <- 1E4

# Spatial random effects
## Simulate GRID adjacency matrix
side <- 10
W_grid <- f_grid_adj(side = side)
n_adj <- nrow(W_grid)
## Generation
set.seed(1312)
alpha_prior <- runif(K)

# Tau
set.seed(617)
tau_prior <- rgamma(K, shape = 2, rate = .2)

# MM matrix
## Parameters
w_ord <- c(.5, .35, .15) # Weight of neighbours orders
ord <- length(w_ord) - 1 # Order of neighbours to include
## Grid
set.seed(873)
H_130_grid <- f_mm_mat(W = W_grid, m = 130, ord = ord, w_ord = w_ord,
                       id_vec = rep(1, nrow(W_grid)))$weight_ext
m <- nrow(H_130_grid)
## Reorder memberships randomly
set.seed(1324)
rand_rows <- c(sample(1:n_adj), sample((n_adj + 1):m))
H_130_grid <- H_130_grid[rand_rows,]
## MM inverse
H_inv <- ginv(H_130_grid)
H_inv_100 <- solve(H_130_grid[1:n_adj, 1:n_adj])
H_inv_70 <- ginv(H_130_grid[1:70, ])
## Other memberships
H_100_grid <- H_130_grid[1:n_adj,]
H_70_grid <- H_130_grid[1:70,]

# Linear term
## Generate parameters
set.seed(432)
gamma_prior <- rnorm(K, sd = .7)
# gamma_prior <- rep(0, K)
beta1_prior <- rnorm(K, sd = .7)
beta2_prior <- rnorm(K, sd = .7)
## Generate min-max covariates
set.seed(34324)
X_raw_grid <- cbind(rnorm(nrow(W_grid)), rnorm(nrow(W_grid)))
X_minmax <- apply(
  X_raw_grid, 2, function(x) (x - min(x))/diff(range(x))
)
## Choose covariates
X_cov <- X_minmax
### Membership level covariates
X_mm <- apply(X_cov, 2, function(x) H_130_grid %*% x)

## Generate offsets
set.seed(3321)
off_130_grid <- rpois(m, lambda = 20)
off_100_grid <- off_130_grid[1:n_adj]
off_70_grid <- off_130_grid[1:70]
#

#### DATA GENERATION - POST ####

# Spatial random effects
phi_grid <- matrix(NA, nrow = K, ncol = nrow(W_grid))
set.seed(223)
for(i in 1:K) phi_grid[i,] <- f_phi_car(
  W = W_grid, alpha = alpha_prior[i], tau = tau_prior[i], B = 1
)

## Global covariate effect
covar_effect <- t(apply(
  cbind(beta1_prior, beta2_prior), 1, function(x) X_cov %*% x
)) + gamma_prior

# Relative risks
## Areal
l_RR_grid <- covar_effect + phi_grid
## Membership
l_RRmm_grid <- apply(l_RR_grid, 1, function(x) H_130_grid %*% x) %>% t()
l_RRmm_100_grid <- l_RRmm_grid[, 1:n_adj]
l_RRmm_70_grid <- l_RRmm_grid[, 1:70]

# Generate outcomes
set.seed(2224)
y_grid_130 <- apply(
  l_RRmm_grid, 1, function(x)  rpois(m, lambda = off_130_grid*exp(x))
) %>% t()
y_grid_100 <- y_grid_130[, 1:n_adj] ; y_grid_70 <- y_grid_130[, 1:70]
#

#### SAVE RESULTS - POST ####

# Data
## 130
d_grid_130 <- list(
  # Phi
  alpha = alpha_prior, tau = tau_prior, phi = phi_grid,
  # Linear term
  X_cov = X_cov, X_mm = X_mm,
  gamma = gamma_prior, beta = cbind(beta1_prior, beta2_prior),
  # Adjacency and membership
  W = W_grid, H = H_130_grid, H_inv = H_inv,
  # Relative risks
  l_RR = l_RR_grid, l_RRmm = l_RRmm_grid,
  # Outcomes
  y = y_grid_130, off_grid = off_130_grid
)
## 100
d_grid_100 <- list(
  # Phi
  alpha = alpha_prior, tau = tau_prior, phi = phi_grid, 
  # Linear term
  X_cov = X_cov, X_mm = X_mm[1:n_adj, ],
  gamma = gamma_prior, beta = cbind(beta1_prior, beta2_prior),
  # Adjacency and membership
  W = W_grid, H = H_100_grid, H_inv = H_inv_100,
  # Relative risks
  l_RR = l_RR_grid, l_RRmm = l_RRmm_100_grid,
  # Outcomes
  y = y_grid_100, off_grid = off_100_grid
)
## 70
d_grid_70 <- list(
  # Phi
  alpha = alpha_prior, tau = tau_prior, phi = phi_grid,
  # Linear term
  X_cov = X_cov, X_mm = X_mm[1:70, ],
  gamma = gamma_prior, beta = cbind(beta1_prior, beta2_prior),
  # Adjacency and membership
  W = W_grid, H = H_70_grid, H_inv = H_inv_70,
  # Relative risks
  l_RR = l_RR_grid, l_RRmm = l_RRmm_70_grid,
  # Outcomes
  y = y_grid_70, off_grid = off_70_grid
)

# Save data and models
save(d_grid_130, d_grid_100, d_grid_70,
     file = "Data/Data_INV_POST_SBC_POST.Rdata")
#

#### DATA GENERATION - INVERSE ####

# Spatial random effects - membership level
## 100
W_t <- t(H_inv_100) %*% W_grid %*% H_inv_100
D_t <- t(H_inv_100) %*% diag(rowSums(W_grid)) %*% H_inv_100
set.seed(3212)
phi_grid_100_mm <- matrix(NA, nrow = K, ncol = n_adj)
for(i in 1:K){
  phi_grid_100_mm[i,] <- mvrnorm(
    n = 1, mu = rep(0, n_adj),
    Sigma = solve(tau_prior[i] * (D_t - alpha_prior[i] * W_t))
  )
}
### Areal level
phi_grid_100 <- apply(phi_grid_100_mm, 1, function(x) H_inv_100 %*% x) %>% t()
## 70
W_t_70 <- t(H_inv_70) %*% W_grid %*% H_inv_70
D_t_70 <- t(H_inv_70) %*% diag(rowSums(W_grid)) %*% H_inv_70
set.seed(3212)
phi_grid_70_mm <- matrix(NA, nrow = K, ncol = 70)
for(i in 1:K){
  phi_grid_70_mm[i,] <- mvrnorm(
    n = 1, mu = rep(0, 70),
    Sigma = solve(tau_prior[i] * (D_t_70 - alpha_prior[i] * W_t_70))
  )
}
### Areal level
phi_grid_70 <- apply(phi_grid_70_mm, 1, function(x) H_inv_70 %*% x) %>% t()


# Covar effect
## 130
covar_effect_130_mm <- t(apply(
  cbind(beta1_prior, beta2_prior), 1, function(x) X_mm %*% x
)) + gamma_prior
## 100
covar_effect_100_mm <- covar_effect_130_mm[, 1:n_adj]
## 70
covar_effect_70_mm <- covar_effect_130_mm[, 1:70]
## 100 areal
covar_effect_100_inv <- t(apply(
  cbind(beta1_prior, beta2_prior), 1, function(x) X_cov %*% x
)) + gamma_prior

# Membership
## 100
l_RRmm_grid_100_inv <- covar_effect_100_mm + phi_grid_100_mm
## 70
l_RRmm_grid_70_inv <- covar_effect_70_mm + phi_grid_70_mm

# Areal
## 100
l_RR_grid_100 <- covar_effect_100_inv + phi_grid_100
## 70
l_RR_grid_70 <- covar_effect_100_inv + phi_grid_70

# Generate outcomes
## 100
set.seed(2224)
y_grid_100_inv <- apply(
  l_RRmm_grid_100_inv, 1, function(x)  rpois(n_adj, lambda = off_100_grid*exp(x))
) %>% t()
## 70
set.seed(2224)
y_grid_70_inv <- apply(
  l_RRmm_grid_70_inv, 1, function(x)  rpois(n_adj, lambda = off_70_grid*exp(x))
) %>% t()
#

#### SAVE RESULTS - INVERSE ####

# True parameters
## 100
d_grid_100 <- list(
  # Phi
  alpha = alpha_prior, tau = tau_prior, phi = phi_grid_100, 
  # Linear term
  X_cov = X_cov, X_mm = X_mm[1:n_adj,],
  gamma = gamma_prior, beta = cbind(beta1_prior, beta2_prior),
  # Adjacency and membership
  W = W_grid, H = H_100_grid, H_inv = H_inv_100,
  # Relative risks
  l_RR = l_RR_grid_100, l_RRmm = l_RRmm_grid_100_inv,
  # Outcomes
  y = y_grid_100_inv, off_grid = off_100_grid
)
## 70
d_grid_70 <- list(
  # Phi
  alpha = alpha_prior, tau = tau_prior, phi = phi_grid_70, 
  # Linear term
  X_cov = X_cov, X_mm = X_mm[1:n_adj,],
  gamma = gamma_prior, beta = cbind(beta1_prior, beta2_prior),
  # Adjacency and membership
  W = W_grid, H = H_70_grid, H_inv = H_inv_70,
  # Relative risks
  l_RR = l_RR_grid_70, l_RRmm = l_RRmm_grid_70_inv,
  # Outcomes
  y = y_grid_70_inv, off_grid = off_70_grid
)

# Save data and models
save(d_grid_100, d_grid_70,
     file = "Data/Data_INV_POST_SBC_INV.Rdata")
#
