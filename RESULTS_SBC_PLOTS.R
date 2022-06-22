#### EXTERNAL PARAMETERS ####

sim_id <- "check"
dim_par <- 100
param <- "post"
# 

#### SETUP ####

# TODO: SET WORKING DIRECTORY
setwd(file.path(path))
source("Functions.R")

# Useful paths
## Images
path_save_img <- file.path("Results/Plots")
## Results
path_res <- file.path("Results/Posterior_samples")
## Control for missing sims
file_list <- list.files(path_res)

# Load data
if(param == "inv"){
  load(file.path(path, "Data/Data_INV_POST_SBC_INV.Rdata"))
} else{
  load(file.path(path, "Data/Data_INV_POST_SBC_POST.Rdata"))
}
## Load single file for parameters (take the first one in the folder)
load(file.path(path_res, file_list[1]))

# Libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(xtable)

## Function to trim posterior rank results
rows_na <- function(mat){
  mat[apply(mat, 1, function(x) sum(is.na(x))) == 0, ]
}
#

#### COMMON PARAMETERS ####

d_grid <- get(paste("d_grid_", dim_par, sep = ""))

# Number of iterations
B <- l_iter$niter*l_iter$nchains/2
par_thinning <- l_iter$par_thinning

## Covariates
X_cov <- d_grid$X_cov
X_mm <- d_grid$X_mm
## Adjacency matrix
W_grid <- d_grid$W
## MM matrix total
H <- d_grid$H
## Total number of simulations
K <- nrow(d_grid$y)

# Global parameters
## General dimensions
n <- nrow(d_grid$W) ; m <- nrow(d_grid$H)
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
#

#### LOAD RESULTS - SBC ####

# Base matrices
par_names <- c(pars_mean, phi_names, RR_names, rr_mm_names)
## INVERSE
sbc_inv <- matrix(NA, nrow = K, ncol = length(par_names) + 1)
colnames(sbc_inv) <- c(par_names, "Rhat")
## POST
sbc_post <- matrix(NA, nrow = K, ncol = length(par_names) + 1)
colnames(sbc_post) <- c(par_names, "Rhat")

for(i in 1:K){
  if( # Check if file is present in the folder
    !paste(
      "sbc_", i, "_", sim_id, "_", dim_par, "_", param, ".RData", sep = ""
    ) %in% file_list
  ){
    next
  }
  
  # Load specific simulation index
  load(file.path(
    path_res, 
    paste("sbc_", i, "_", sim_id, "_", dim_par, "_", param, ".RData", sep = "")
  ))
  
  ## INVERSE
  if(dim_par != 130){
    sbc_inv[i, par_names] <- results_inv[par_names, "rank_stat"]
    sbc_inv[i, "Rhat"] <- max(round(results_inv[, "Rhat"], 4), na.rm = T)
  }
  ## POST
  sbc_post[i, par_names] <- results_post[par_names, "rank_stat"]
  sbc_post[i, "Rhat"] <- max(round(results_post[, "Rhat"], 4), na.rm = T)
}

# Trim - REMOVE UNCONVERGED SIMULATIONS BY LOOKING AT RHATS ABOVE 1.01
Rhat <- 1.01
## POST
unconv_post <- which(sbc_post[,"Rhat"] > Rhat)
write(
  unconv_post,
  file = file.path(
    "Results", paste("unconv_", dim_par,"_", param, "_post.txt", sep = "")
  )
)
## INVERSE
unconv_inv <- which(sbc_inv[,"Rhat"] > Rhat)
write(
  unconv_inv,
  file = file.path(
    "Results", paste("unconv_", dim_par,"_", param, "_inv.txt", sep = "")
  )
)
## NA - remove NA rows 
if(dim_par != 130) sbc_inv <- rows_na(sbc_inv)
sbc_post <- rows_na(sbc_post)
## Rhat - REMOVE UNCONVERGED MCMCs
sbc_inv <- sbc_inv[which(sbc_inv[,"Rhat"] <= Rhat),]
sbc_post <- sbc_post[which(sbc_post[,"Rhat"] <= Rhat),]
#

#### COMMON FOR PLOTS ####

## Interval
global_p_cov <- .025
## Number of thinned samples
thinned_B <- length(seq(1, B, by = par_thinning))
## Colours
p_colours <- data.frame(x = 1:128, y = 0, col = 1:128) %>% 
  ggplot(aes(x = x, y = y, colour = as.character(col))) + geom_point(size = 2)
v_colours <- ggplot_build(p_colours)$data[[1]]$colour
## Titles
quant <- c("y", "rr", "phi")
## title lookup for plot titles
title_chr <- cbind(
  index = c("cal", "sbc", "inv", "post", "post_un", quant,
            "alpha", "tau", "gamma"),
  title = c("Calibration", "SBC", "Inverse", "Post", "Post (Unstructured)",
            "Observed (y)", "rho", "phi",
            "alpha", "tau", "gamma")
)

# Common arguments
args_common <- list(p_cov = global_p_cov, gg_base = base_plot)
#

#### PRINT PLOTS - SBC MANUAL - MCMC: INV ####

# Fix indeces
p_type <- c("sbc")
quant_sbc_ind <- c(
  "rr", "phi", "gamma", "tau", "alpha"
)
quant_sbc_names <- list(
  rr = RR_names,
  phi = phi_names,
  rrmm = rr_mm_names,
  alpha = "alpha", tau = "tau", gamma = "gamma", 
  `beta[1]` = "beta[1]", `beta[2]` = "beta[2]"
)

# Allowed MCMC parameterisations
if(dim_par != 130){
  samp_par <- c("inv", "post")
} else{
  samp_par <- "post"
}

## Scaling functions
f_scale <- function(x) x
f_inv_scale <- function(x) x

# INVERSE SAMPLING
if(dim_par != 130){
  ## RR MM
  current_plot <- do.call(
    get(paste("f_", p_type, "_plot", sep = "")),
    c(
      args_common, 
      list(
        d_rank = sbc_inv,
        L_thin = thinned_B, vec_colours = v_colours,
        v_pars = quant_sbc_names[["rrmm"]]
      )
    )
  )
  ggsave(
    file = file.path(
      path_save_img,
      paste("sbc_inv_rrmm_", dim_par, ".png", sep = "")
    ),
    plot = current_plot + 
      ggtitle(
        bquote(
          "Data generation: " ~ 
            .(title_chr[which(title_chr[, "index"] == param), "title"]) ~ 
            " -  MCMC: " ~
            .(title_chr[which(title_chr[, "index"] == "inv"), "title"]) ~
            " - " ~ tilde(rho) ~ " - " ~
            " (m = " ~  .(dim_par) ~")"
          
        )
      ) + theme(
        axis.text.x = element_text(face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 15),
        plot.title = element_text(size = 15, face = "bold"),
      ),
    device = "png", dpi = "retina"
  )
  
  ## beta 1
  current_plot <- do.call(
    get(paste("f_", p_type, "_plot", sep = "")),
    c(
      args_common, 
      list(
        d_rank = sbc_inv,
        L_thin = thinned_B, vec_colours = v_colours,
        v_pars = quant_sbc_names[["beta[1]"]]
      )
    )
  )
  ggsave(
    file = file.path(
      path_save_img,
      paste("sbc_inv_beta[1]_", dim_par, ".png", sep = "")
    ),
    plot = current_plot + 
      ggtitle(
        bquote(
          "Data generation: " ~ 
            .(title_chr[which(title_chr[, "index"] == param), "title"]) ~ 
            " -  MCMC: " ~
            .(title_chr[which(title_chr[, "index"] == "inv"), "title"]) ~
            " - " ~ beta[1] ~ " - " ~
            " (m = " ~  .(dim_par) ~")"
          
        )
      ) + theme(
        axis.text.x = element_text(face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 15),
        plot.title = element_text(size = 15, face = "bold"),
      ),
    device = "png", dpi = "retina"
  )
  
  ## beta 2
  current_plot <- do.call(
    get(paste("f_", p_type, "_plot", sep = "")),
    c(
      args_common, 
      list(
        d_rank = sbc_inv,
        L_thin = thinned_B, vec_colours = v_colours,
        v_pars = quant_sbc_names[["beta[2]"]]
      )
    )
  )
  ggsave(
    file = file.path(
      path_save_img,
      paste("sbc_inv_beta[2]_", dim_par, ".png", sep = "")
    ),
    plot = current_plot + 
      ggtitle(
        bquote(
          "Data generation: " ~ 
            .(title_chr[which(title_chr[, "index"] == param), "title"]) ~ 
            " -  MCMC: " ~
            .(title_chr[which(title_chr[, "index"] == "inv"), "title"]) ~
            " - " ~ beta[2] ~ " - " ~
            " (m = " ~  .(dim_par) ~")"
          
        )
      ) + theme(
        axis.text.x = element_text(face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 15),
        plot.title = element_text(size = 15, face = "bold"),
      ),
    device = "png", dpi = "retina"
  )
}
#

#### PRINT PLOTS- SBC MANUAL - MCMC: POST ####

current_plot <- do.call(
  get(paste("f_", p_type, "_plot", sep = "")),
  c(
    args_common, 
    list(
      d_rank = sbc_post,
      L_thin = thinned_B, vec_colours = v_colours,
        v_pars = quant_sbc_names[["rrmm"]]
    )
  )
)
ggsave(
  file = file.path(
    path_save_img,
    paste("sbc_post_rrmm_", dim_par, ".png", sep = "")
  ),
  plot = current_plot + 
    ggtitle(
      bquote(
        "Data generation: " ~ 
          .(title_chr[which(title_chr[, "index"] == param), "title"]) ~ 
          " -  MCMC: " ~
          .(title_chr[which(title_chr[, "index"] == "post"), "title"]) ~
          " - " ~ tilde(rho) ~ " - " ~
          " (m = " ~  .(dim_par) ~")"
        
      )
    ) + theme(
      axis.text.x = element_text(face = "bold", size = 15),
      axis.text.y = element_text(face = "bold", size = 15),
      plot.title = element_text(size = 15, face = "bold"),
    ),
  device = "png", dpi = "retina"
)

## beta 1
current_plot <- do.call(
  get(paste("f_", p_type, "_plot", sep = "")),
  c(
    args_common, 
    list(
      d_rank = sbc_post,
      L_thin = thinned_B, vec_colours = v_colours,
      v_pars = quant_sbc_names[["beta[1]"]]
    )
  )
)
ggsave(
  file = file.path(
    path_save_img,
    paste("sbc_post_beta[1]_", dim_par, ".png", sep = "")
  ),
  plot = current_plot + 
    ggtitle(
      bquote(
        "Data generation: " ~ 
          .(title_chr[which(title_chr[, "index"] == param), "title"]) ~ 
          " -  MCMC: " ~
          .(title_chr[which(title_chr[, "index"] == "post"), "title"]) ~
          " - " ~ beta[1] ~ " - " ~
          " (m = " ~  .(dim_par) ~")"
        
      )
    ) + theme(
      axis.text.x = element_text(face = "bold", size = 15),
      axis.text.y = element_text(face = "bold", size = 15),
      plot.title = element_text(size = 15, face = "bold"),
    ),
  device = "png", dpi = "retina"
)

## beta 2
current_plot <- do.call(
  get(paste("f_", p_type, "_plot", sep = "")),
  c(
    args_common, 
    list(
      d_rank = sbc_post,
      L_thin = thinned_B, vec_colours = v_colours,
      v_pars = quant_sbc_names[["beta[2]"]]
    )
  )
)
ggsave(
  file = file.path(
    path_save_img,
    paste("sbc_post_beta[2]_", dim_par, ".png", sep = "")
  ),
  plot = current_plot + 
    ggtitle(
      bquote(
        "Data generation: " ~ 
          .(title_chr[which(title_chr[, "index"] == param), "title"]) ~ 
          " -  MCMC: " ~
          .(title_chr[which(title_chr[, "index"] == "post"), "title"]) ~
          " - " ~ beta[2] ~ " - " ~
          " (m = " ~  .(dim_par) ~")"
        
      )
    ) + theme(
      axis.text.x = element_text(face = "bold", size = 15),
      axis.text.y = element_text(face = "bold", size = 15),
      plot.title = element_text(size = 15, face = "bold"),
    ),
  device = "png", dpi = "retina"
)
#

#### BIAS TOTAL QUANTITES ####

# Create total quantities - POST
## MM
base_mat <- matrix(NA, nrow = K, ncol = m)
colnames(base_mat) <- rr_mm_names
bias_mm_post <- base_mat ; abs_mm_post <- base_mat ; rmse_mm_post <- base_mat
## RR
base_mat <- matrix(NA, nrow = K, ncol = n)
colnames(base_mat) <- RR_names
bias_rr_post <- base_mat ; abs_rr_post <- base_mat ; rmse_rr_post <- base_mat
## PHI
base_mat <- matrix(NA, nrow = K, ncol = n)
colnames(base_mat) <- phi_names
bias_phi_post <- base_mat ; abs_phi_post <- base_mat ; rmse_phi_post <- base_mat
## Mean
base_mat <- matrix(NA, nrow = K, ncol = length(pars_mean))
colnames(base_mat) <- pars_mean
bias_pars_post <- base_mat ; abs_pars_post <- base_mat ; rmse_pars_post <- base_mat

# Create total quantities - INV
if(dim_par != 130){
  ## MM
  base_mat <- matrix(NA, nrow = K, ncol = m)
  colnames(base_mat) <- rr_mm_names
  bias_mm_inv <- base_mat ; abs_mm_inv <- base_mat ; rmse_mm_inv <- base_mat
  ## RR
  base_mat <- matrix(NA, nrow = K, ncol = n)
  colnames(base_mat) <- RR_names
  bias_rr_inv <- base_mat ; abs_rr_inv <- base_mat ; rmse_rr_inv <- base_mat
  ## PHI
  base_mat <- matrix(NA, nrow = K, ncol = n)
  colnames(base_mat) <- phi_names
  bias_phi_inv <- base_mat ; abs_phi_inv <- base_mat ; rmse_phi_inv <- base_mat
  ## Mean
  base_mat <- matrix(NA, nrow = K, ncol = length(pars_mean))
  colnames(base_mat) <- pars_mean
  bias_pars_inv <- base_mat ; abs_pars_inv <- base_mat ; rmse_pars_inv <- base_mat
}
#

#### BIAS - RETRIEVE RESULTS ####

for(i in 1:K){
  if(
    !paste("sbc_", i, "_", sim_id, "_", dim_par, "_", param, ".RData", sep = "") %in% file_list
  ){
    next
  }
  
  load(file.path(
    path_res, paste("sbc_", i, "_", sim_id, "_", dim_par, "_", param, ".RData", sep = "")
  ))
  
  ## INVERSE
  if(dim_par != 130){
    ### BIAS
    bias_mm_inv[i, rr_mm_names] <- results_inv[rr_mm_names, "bias"]
    bias_rr_inv[i, RR_names] <- results_inv[RR_names, "bias"]
    bias_phi_inv[i, phi_names] <- results_inv[phi_names, "bias"]
    bias_pars_inv[i, pars_mean] <- results_inv[pars_mean, "bias"]
    ### ABS
    abs_mm_inv[i, rr_mm_names] <- results_inv[rr_mm_names, "abs_bias"]
    abs_rr_inv[i, RR_names] <- results_inv[RR_names, "abs_bias"]
    abs_phi_inv[i, phi_names] <- results_inv[phi_names, "abs_bias"]
    abs_pars_inv[i, pars_mean] <- results_inv[pars_mean, "abs_bias"]
    ### RMSE
    rmse_mm_inv[i, rr_mm_names] <- results_inv[rr_mm_names, "RMSE"]
    rmse_rr_inv[i, RR_names] <- results_inv[RR_names, "RMSE"]
    rmse_phi_inv[i, phi_names] <- results_inv[phi_names, "RMSE"]
    rmse_pars_inv[i, pars_mean] <- results_inv[pars_mean, "RMSE"]
  }
  ## POST
  ### BIAS
  bias_mm_post[i, rr_mm_names] <- results_post[rr_mm_names, "bias"]
  bias_rr_post[i, RR_names] <- results_post[RR_names, "bias"]
  bias_phi_post[i, phi_names] <- results_post[phi_names, "bias"]
  bias_pars_post[i, pars_mean] <- results_post[pars_mean, "bias"]
  ### ABS
  abs_mm_post[i, rr_mm_names] <- results_post[rr_mm_names, "abs_bias"]
  abs_rr_post[i, RR_names] <- results_post[RR_names, "abs_bias"]
  abs_phi_post[i, phi_names] <- results_post[phi_names, "abs_bias"]
  abs_pars_post[i, pars_mean] <- results_post[pars_mean, "abs_bias"]
  ### RMSE
  rmse_mm_post[i, rr_mm_names] <- results_post[rr_mm_names, "RMSE"]
  rmse_rr_post[i, RR_names] <- results_post[RR_names, "RMSE"]
  rmse_phi_post[i, phi_names] <- results_post[phi_names, "RMSE"]
  rmse_pars_post[i, pars_mean] <- results_post[pars_mean, "RMSE"]
  
}

# Trim
if(dim_par != 130){
  ## BIAS
  bias_mm_inv <- rows_na(bias_mm_inv)
  bias_rr_inv <- rows_na(bias_rr_inv)
  bias_phi_inv <- rows_na(bias_phi_inv)
  bias_pars_inv <- rows_na(bias_pars_inv)
  ## ABS
  abs_mm_inv <- rows_na(abs_mm_inv)
  abs_rr_inv <- rows_na(abs_rr_inv)
  abs_phi_inv <- rows_na(abs_phi_inv)
  abs_pars_inv <- rows_na(abs_pars_inv)
  ## RMSE
  rmse_mm_inv <- rows_na(rmse_mm_inv)
  rmse_rr_inv <- rows_na(rmse_rr_inv)
  rmse_phi_inv <- rows_na(rmse_phi_inv)
  rmse_pars_inv <- rows_na(rmse_pars_inv)
}
## BIAS
bias_mm_post <- rows_na(bias_mm_post)
bias_rr_post <- rows_na(bias_rr_post)
bias_phi_post <- rows_na(bias_phi_post)
bias_pars_post <- rows_na(bias_pars_post)
## ABS
abs_mm_post <- rows_na(abs_mm_post)
abs_rr_post <- rows_na(abs_rr_post)
abs_phi_post <- rows_na(abs_phi_post)
abs_pars_post <- rows_na(abs_pars_post)
## RMSE
rmse_mm_post <- rows_na(rmse_mm_post)
rmse_rr_post <- rows_na(rmse_rr_post)
rmse_phi_post <- rows_na(rmse_phi_post)
rmse_pars_post <- rows_na(rmse_pars_post)
#

#### CREATE TABLES SUMMARIES - POST ####

# BIAS
## CAR
bias_pars_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) -> 
  tab_bias_pars
## RR
bias_rr_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) -> 
  tab_bias_rr_post
## MEMBERSHIP
bias_mm_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) ->
  tab_bias_mm_post
## PHI
bias_phi_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) ->
  tab_bias_phi_post

# ABSOLUTE
## CAR
abs_pars_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) -> 
  tab_abs_pars
## RR
abs_rr_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) -> 
  tab_abs_rr_post
## MEMBERSHIP
abs_mm_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) ->
  tab_abs_mm_post
## PHI
abs_phi_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) ->
  tab_abs_phi_post

# RMSE
## CAR
rmse_pars_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) -> 
  tab_rmse_pars
## RR
rmse_rr_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) -> 
  tab_rmse_rr_post
## MEMBERSHIP
rmse_mm_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) ->
  tab_rmse_mm_post
## PHI
rmse_phi_post %>% as.data.frame() %>% 
  pivot_longer(
    cols = everything(), names_to = "names", values_to = "values"
  ) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), median = median(values),
            ql = quantile(values, .025), qh = quantile(values, .975)) ->
  tab_rmse_phi_post
#

#### CREATE TABLES SUMMARIES - INV ####

if(dim_par != 130){
  # BIAS
  ## CAR
  bias_pars_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) -> 
    tab_bias_pars
  ## RR
  bias_rr_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) -> 
    tab_bias_rr_inv
  ## MEMBERSHIP
  bias_mm_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) ->
    tab_bias_mm_inv
  ## PHI
  bias_phi_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) ->
    tab_bias_phi_inv
  
  # ABSOLUTE
  ## CAR
  abs_pars_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) -> 
    tab_abs_pars
  ## RR
  abs_rr_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) -> 
    tab_abs_rr_inv
  ## MEMBERSHIP
  abs_mm_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) ->
    tab_abs_mm_inv
  ## PHI
  abs_phi_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) ->
    tab_abs_phi_inv
  
  # RMSE
  ## CAR
  rmse_pars_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) -> 
    tab_rmse_pars
  ## RR
  rmse_rr_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) -> 
    tab_rmse_rr_inv
  ## MEMBERSHIP
  rmse_mm_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) ->
    tab_rmse_mm_inv
  ## PHI
  rmse_phi_inv %>% as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "names", values_to = "values"
    ) %>% 
    group_by(names) %>% 
    summarise(mean = mean(values), median = median(values),
              ql = quantile(values, .025), qh = quantile(values, .975)) ->
    tab_rmse_phi_inv
}
#

#### TABLES - BUILD XTABLES ####

library(xtable)

# Common parameters
## Row names 
r_names_relative <- c(
  "\\multirow{8}{*}{mean}", "", "", "", "", "", "", "",
  "\\multirow{8}{*}{median}", "", "", "", "", "", "", "",
  "\\multirow{8}{*}{$2.5\\%$}", "", "", "", "", "", "", "",
  "\\multirow{8}{*}{$97.5\\%$}", "", "", "", "", "", "", ""
)
par_names <- c(
  "$\\alpha$", "$\\tau$", "$\\gamma$", "$\\beta_1$", "$\\beta_2$",
  "$\\bm{\\phi}$", "$\\bm{\\rho}$", "$\\tilde{\\bm{\\rho}}$"
)
## Length of summary statistics considered
l_summary <- 4

# BIAS
## Parameters of interest
### alpha tau gamma beta_1 beta_2 phi rho rho_tilde
## Create base table
tab_bias_tot <- matrix(NA, ncol = 9, nrow = 8*l_summary)
colnames(tab_bias_tot) <- c(
  "$70$", "$100$", "$130$", "$70$", "$100$", "$130$"
)
## Add data
### CAR
tab_bias_tot[seq(1, 8*l_summary, by = 8),] <- 
  tab_bias_pars %>% filter(names == "alpha") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_bias_tot[seq(2, 8*l_summary, by = 8),] <- 
  tab_bias_pars %>% filter(names == "tau") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_bias_tot[seq(3, 8*l_summary, by = 8),] <- 
  tab_bias_pars %>% filter(names == "gamma") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_bias_tot[seq(4, 8*l_summary, by = 8),] <- 
  tab_bias_pars %>% filter(names == "beta[1]") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_bias_tot[seq(5, 8*l_summary, by = 8),] <- 
  tab_bias_pars %>% filter(names == "beta[2]") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
### PHI
tab_bias_tot[seq(6, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_bias_phi[,-1]))
### RR
tab_bias_tot[seq(7, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_bias_rr[,-1]))
### MM
tab_bias_tot[seq(8, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_bias_mm[,-1]))
tab_bias_tot <- as.data.frame(tab_bias_tot)
## xtable
print(xtable(cbind(rep(par_names, 4), tab_bias_tot)), 
      booktabs = T,
      # sanitize.rownames.function = function(x) paste(r_names_relative),
      include.rownames=F,
      sanitize.colnames.function = function(x) 
        paste(c("", colnames(tab_bias_tot))),
      sanitize.text.function=identity
)

# ABSOLUTE
## Parameters of interest
### alpha tau gamma beta_1 beta_2 phi rho rho_tilde
## Create base table
tab_abs_tot <- matrix(NA, ncol = 9, nrow = 8*l_summary)
colnames(tab_abs_tot) <- c("$48$", "$54$", "$61$", "$64$", "$77$", "$90$",
                           "$102$", "$115$", "$128$")
## Add data
### CAR
tab_abs_tot[seq(1, 8*l_summary, by = 8),] <- 
  tab_abs_pars %>% filter(names == "alpha") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_abs_tot[seq(2, 8*l_summary, by = 8),] <- 
  tab_abs_pars %>% filter(names == "tau") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_abs_tot[seq(3, 8*l_summary, by = 8),] <- 
  tab_abs_pars %>% filter(names == "gamma") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_abs_tot[seq(4, 8*l_summary, by = 8),] <- 
  tab_abs_pars %>% filter(names == "beta[1]") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_abs_tot[seq(5, 8*l_summary, by = 8),] <- 
  tab_abs_pars %>% filter(names == "beta[2]") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
### PHI
tab_abs_tot[seq(6, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_abs_phi[,-1]))
### RR
tab_abs_tot[seq(7, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_abs_rr[,-1]))
### MM
tab_abs_tot[seq(8, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_abs_mm[,-1]))
tab_abs_tot <- as.data.frame(tab_abs_tot)
## xtable
print(xtable(cbind(rep(par_names, 4), tab_abs_tot)), 
      booktabs = T,
      # sanitize.rownames.function = function(x) paste(r_names_relative),
      include.rownames=F,
      sanitize.colnames.function = function(x) 
        paste(c("", colnames(tab_abs_tot))),
      sanitize.text.function=identity
)

# RMSE
## Parameters of interest
### alpha tau gamma beta_1 beta_2 phi rho rho_tilde
## Create base table
tab_rmse_tot <- matrix(NA, ncol = 9, nrow = 8*l_summary)
colnames(tab_rmse_tot) <- c("$48$", "$54$", "$61$", "$64$", "$77$", "$90$",
                            "$102$", "$115$", "$128$")
## Add data
### CAR
tab_rmse_tot[seq(1, 8*l_summary, by = 8),] <- 
  tab_rmse_pars %>% filter(names == "alpha") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_rmse_tot[seq(2, 8*l_summary, by = 8),] <- 
  tab_rmse_pars %>% filter(names == "tau") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_rmse_tot[seq(3, 8*l_summary, by = 8),] <- 
  tab_rmse_pars %>% filter(names == "gamma") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_rmse_tot[seq(4, 8*l_summary, by = 8),] <- 
  tab_rmse_pars %>% filter(names == "beta[1]") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
tab_rmse_tot[seq(5, 8*l_summary, by = 8),] <- 
  tab_rmse_pars %>% filter(names == "beta[2]") %>% ungroup() %>% 
  select(-c(memb, names)) %>% t() %>% as.matrix()
### PHI
tab_rmse_tot[seq(6, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_rmse_phi[,-1]))
### RR
tab_rmse_tot[seq(7, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_rmse_rr[,-1]))
### MM
tab_rmse_tot[seq(8, 8*l_summary, by = 8),] <- 
  as.matrix(t(tab_rmse_mm[,-1]))
tab_rmse_tot <- as.data.frame(tab_rmse_tot)
## xtable
print(xtable(cbind(rep(par_names, 4), tab_rmse_tot)), 
      booktabs = T,
      # sanitize.rownames.function = function(x) paste(r_names_relative),
      include.rownames=F,
      sanitize.colnames.function = function(x) 
        paste(c("", colnames(tab_rmse_tot))),
      sanitize.text.function=identity
)
#




