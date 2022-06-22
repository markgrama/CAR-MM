#### GENERAL SETUP ####

# TODO: SET WORKING DIRECTORY
setwd(path)

# Load Data
load("DataAnalysis/SEL_data_DEF.Rdata")

# Load package
## Plot
# library(plot.matrix)
# library(viridis)
## Data handling
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(stringr)
# library(forcats)
## Sampling
# library(MASS)
## Correlation
# library(corpcor)
# library(correlation)
# library(MASS)
## Spatial
# library(rgeos)
# library(rgdal)
## Matrices
# library(Matrix)
## Stan
library(rstan)
# library(shinystan)
library(loo)
rstan_options(auto_write = TRUE)      
options(mc.cores = parallel::detectCores())
## MCMC
# library(bayesplot)

# Load functions
# source("Scripts/CARMM_Functions.R")
# source("Scripts/Posterior_checks.R")
# source("Scripts/CARMM_Structure_Functions.R")

# Posterior fits paths
res_path <- file.path("DataAnalysis")

# Data
## Centre covariates
X_cent <- apply(
  X_msoa[, c("perc_sa", "imd")], 2, function(x) (x - min(x))/diff(range(x))
)
## Bind
d_sel <- list(
  # Number of areas
  n = nrow(W_sel),
  # Covariates
  k = ncol(X_cent),
  X_cov = X_cent,
  # Adjacecncy
  W_n = sum(W_sel) / 2,
  # number of neighbor pairs
  W = W_sel,
  # Multiple membership
  m = nrow(H_sel),
  H = H_sel,
  # Outcome
  y = outcomes$prev,
  log_offset = log(outcomes$exp_prev)
)
#

#### HMC PARAMTERS ####

niter <- 1E4
nchains <- 5
#

#### GLM-MM ####

mod_glm <- stan_model(
  file = "DataAnalysis/Models/GLMM_COV_NB.stan"
)

time_glm <- Sys.time()
fit_glm <- sampling(
  mod_glm, data = d_sel, iter = niter, chains = nchains,
  control = list(adapt_delta = .99, max_treedepth = 15)
)
time_glm <- Sys.time() - time_glm
## Results
pars_vec <- c('gamma', 'beta','psi', 'ppp', 'lp__')
print(fit_glm, pars = pars_vec, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))

# Save
save(
  time_glm, fit_glm,
  file = paste(
    "Results/Fit_GLM-MM_", Sys.Date(), ".RData", sep = ""
  )
)
#

#### CAR-MM ####

mod_car <- stan_model(
  file = "Models/CARMM_COV_NB.stan"
)

time_car <- Sys.time()
fit_car <- sampling(
  mod_car, data = d_sel, iter = niter, chains = nchains,
  control = list(adapt_delta = .99, max_treedepth = 15)
)
time_car <- Sys.time() - time_car
## Results
pars_vec <- c('gamma', 'beta','psi', 'ppp', 'lp__')
print(fit_car, pars = pars_vec, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))

# Save
save(
  time_car, fit_car,
  file = paste(
    "DataAnalysis/Fit_CAR-MM_", Sys.Date(), ".RData", sep = ""
  )
)
#

#### ICAR-MM ####

mod_icar <- stan_model(
  file = "DataAnalysis/Models/ICARMM_COV_NB.stan"
)

time_icar <- Sys.time()
fit_icar <- sampling(
  mod_icar, data = d_sel, iter = niter, chains = nchains,
  control = list(adapt_delta = .99, max_treedepth = 15)
)
time_icar <- Sys.time() - time_icar
## Results
pars_vec <- c('gamma', 'beta','psi', 'ppp', 'lp__')
print(fit_icar, pars = pars_vec, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))

# Save
save(
  time_icar, fit_icar,
  file = paste(
    "Results/Fit_ICAR-MM_", Sys.Date(), ".RData", sep = ""
  )
)
#
  
#### LOO ####

# GLM
## Likelihood
lik_glm <- extract_log_lik(
  fit_glm, parameter_name = "log_lik", merge_chains = FALSE
)
## Reff
r_eff_glm <- relative_eff(exp(lik_glm), cores = 2)
## LOO
loo_glm <- loo(lik_glm, r_eff = r_eff_glm, cores = 2)

# GLM - UNSTRUCTURED
## Likelihood
lik_glm_unst <- extract_log_lik(
  fit_glm_unst, parameter_name = "log_lik", merge_chains = FALSE
)
## Reff
r_eff_glm_unst <- relative_eff(exp(lik_glm_unst), cores = 2)
## LOO
loo_glm_unst <- loo(lik_glm_unst, r_eff = r_eff_glm_unst, cores = 2)

# CAR 
## Likelihood
lik_car <- extract_log_lik(
  fit_car, parameter_name = "log_lik", merge_chains = FALSE
)
## Reff
r_eff_car <- relative_eff(exp(lik_car), cores = 2)
## LOO
loo_car <- loo(lik_car, r_eff = r_eff_car, cores = 2)

# ICAR
## Likelihood
lik_icar <- extract_log_lik(
  fit_icar, parameter_name = "log_lik", merge_chains = FALSE
)
## Reff
r_eff_icar <- relative_eff(exp(lik_icar), cores = 2)
## LOO
loo_icar <- loo(lik_icar, r_eff = r_eff_icar, cores = 2)

# Compare
loo_compare(
  list(
    glm = loo_glm, unst = loo_glm_unst,
    car = loo_car, icar = loo_icar
  )
)
#

#### MIXED REPLICATES ####

# GLM
rep_glm <- f_mix_carmm(
  glm = T, nb = T,
  beta = as.matrix(fit_glm)[,c("beta[1]", "beta[2]")],
  gamma = as.matrix(fit_glm)[, "gamma"],
  covars = d_sel$X_cov, W = d_sel$W, H = d_sel$H,
  psi = as.matrix(fit_glm)[, "psi"],
  off = exp(d_sel$log_offset)
)

# CAR
rep_car <- f_mix_carmm(
  nb = T,
  beta = as.matrix(fit_car)[,c("beta[1]", "beta[2]")],
  gamma = as.matrix(fit_car)[, "gamma"],
  tau = as.matrix(fit_car)[, "tau"],
  alpha = as.matrix(fit_car)[, "alpha"],
  covars = d_sel$X_cov, W = d_sel$W, H = d_sel$H,
  psi = as.matrix(fit_car)[, "psi"],
  off = exp(d_sel$log_offset)
)

# ICAR
rep_icar <- f_mix_carmm(
  nb = T, icar = T,
  beta = as.matrix(fit_icar)[,c("beta[1]", "beta[2]")], 
  gamma = as.matrix(fit_icar)[, "gamma"],
  tau = as.matrix(fit_icar)[, "tau"],
  covars = d_sel$X_cov, W = d_sel$W, H = d_sel$H,
  psi = as.matrix(fit_icar)[, "psi"],
  off = exp(d_sel$log_offset)
)
#

#### SCORING RULES ####

# GLM
rps(y_rep = rep_glm$y_rep, y_obs = d_sel$y)
dss(y_rep = rep_glm$y_rep, y_obs = d_sel$y)

# GLM - UNST
rps(y_rep = rep_glm_unst$y_rep, y_obs = d_sel$y)
dss(y_rep = rep_glm_unst$y_rep, y_obs = d_sel$y)

# CAR
rps(y_rep = rep_car$y_rep, y_obs = d_sel$y)
dss(y_rep = rep_car$y_rep, y_obs = d_sel$y)

# ICAR
rps(y_rep = rep_icar$y_rep, y_obs = d_sel$y)
dss(y_rep = rep_icar$y_rep, y_obs = d_sel$y)
#

#### MAP RR ####

# Names vectors
l_RR_names <- sapply(1:d_sel$n, function(x) paste(
  "l_RR[", x, "]", sep = ""
))

# Posterior Means RR for models
## GLM
rr_sel_glm <- apply(
  as.matrix(fit_glm)[, l_RR_names], 2, 
  function(x) mean(exp(x))
)
## CAR
rr_sel_car <- apply(
  as.matrix(fit_car)[, l_RR_names], 2, 
  function(x) mean(exp(x))
)
## ICAR
rr_sel_icar <- apply(
  as.matrix(fit_icar)[, l_RR_names], 2, 
  function(x) mean(exp(x))
)
## Merge
rr_merge <- data.frame(
  id = rownames(d_sel$W),
  glm = rr_sel_glm, car = rr_sel_car, icar = rr_sel_icar
)
## Complete data merge
sel_shp_df <- fortify(sel_df_shp, region = "msoa11cd")
sel_df <- merge(sel_shp_df %>% filter(id %in% rownames(d_sel$W)),
                rr_merge, by = "id")

# Map - BASE
map_base <- ggplot() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), aspect.ratio = 9/16,
        panel.background = element_blank())
## GLM
map_glm <- map_base +
  geom_polygon(data = sel_df, aes(x = long, y = lat, group = id, fill = glm),
               color = 'black', size = .2) +
  scale_fill_viridis(name = bquote("Relative risk" ~ rho)) + 
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15, face = "bold")
  )
## ICAR
map_icar <- map_base +
  geom_polygon(data = sel_df, aes(x = long, y = lat, group = id, fill = icar),
               color = 'black', size = .2) +
  scale_fill_viridis(name = bquote("Relative risk" ~ rho)) + 
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15, face = "bold")
  )
## CARMM
map_car <- map_base +
  geom_polygon(data = sel_df, aes(x = long, y = lat, group = id, fill = car),
               color = 'black', size = .2) +
  scale_fill_viridis(name = bquote("Relative risk" ~ rho)) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15, face = "bold")
  )
#

#### CALIBRATION - RR ####

# Common quantities
H <- d_sel$H
y_obs <- outcomes$prev
m <- nrow(d_sel$H)

# Mixed PPP
order_mm_stat <- 1:m/(m + 1)
## GLM-MM
ppp_glm <- f_ppp_cont(y = y_obs, yrep = t(rep_glm$y_rep)) 
## CAR-MM 
ppp_car <- f_ppp_cont(y = y_obs, yrep = t(rep_car$y_rep)) 
## ICAR-MM
ppp_icar <- f_ppp_cont(y = y_obs, yrep = t(rep_icar$y_rep)) 

# Plot
cal_nb_plot <- data.frame(
  ord = order_mm_stat, 
  glm = sort(ppp_glm), car = sort(ppp_car), icar = sort(ppp_icar)
) %>% 
  pivot_longer(!ord, names_to = "names", values_to = "values") %>% 
  ggplot(aes(x = ord, y = values, colour = names)) + geom_point() +
  scale_colour_viridis(
    discrete = TRUE, name = "Model", 
    labels = c("GLM-MM", "CAR-MM", "ICAR-MM")
  ) + 
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  ) +
  xlab("Ordered Mixed predictive p-values") +
  ylab("Ordinates of Uniform distribution")

plot(order_mm_stat, sort(ppp_icar_diff))
points(order_mm_stat, sort(ppp_glm_nb), col = "red")
points(order_mm_stat, sort(ppp_icar_mor), col = "blue")
abline(a = 0, b = 1)
#

