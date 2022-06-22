#### LIBRARIES ####

library(ggplot2)
#

#### ADJACENCY MATRIX FOR GRID ####

f_grid_adj <- function(side){
  x.easting <- 1:side
  x.northing <- 1:side
  Grid <- expand.grid(x.easting, x.northing)
  K <- nrow(Grid)
  ## Neighbours matrix: based on sharing a common border
  ### Use manhattan distance to also compute 2nd and 3rd order neighbours
  distance <- as.matrix(dist(Grid, method = "manhattan"))
  W <- array(0, c(K,K))
  W[distance == 1] <- 1
  
  return(W)
}
#

#### MULTIPLE MEMBERSHIP MATRIX SIMULATION ####

f_mm_high_ord <- function(
    W, ord, mat_ord_weight, ind_range
){
  n <- nrow(W)
  # Create base matrix of weights
  weight_ext <- matrix(0, ncol = n, nrow = length(ind_range))
  
  # Compute necessary multiple order of neighbours matrices
  W_cml <- W
  for(i in 3:(ord+1)){ # Starts from 3 bc: 1=identical order ; 2=first order neigh
    # Name for i-th order neighbour matrix
    # name_ind <- paste("W_", i, sep = "")
    
    # Compute all areas reachable in EXACTLY i-1 steps
    W_cml_current <- W %^% (i-1)
    ## Transform in adj matrix
    W_cml_current[W_cml_current != 0] <- 1
    ## Remove unwanted diagonal
    diag(W_cml_current) <- rep(0, nrow(W_cml_current))
    
    # Find exactly i-th order neighbour
    W_now <- W_cml_current - W_cml
    ## Adjust matrix
    W_now[W_now < 0] <- 0
    W_now[W_now != 0] <- 1
    
    # Simulate weights for current neighbour order
    ## Simulate random weights
    W_now[which(W_now != 0)] <- runif(length(which(W_now != 0)))
    for(j in ind_range){
      W_now[j, ] <- (W_now[j,]/sum(W_now[j,]))*mat_ord_weight[j,i]
    }
    
    # Save exactly i-th order neighbour
    # assign(name_ind, W_now)W
    
    # Save global matrix
    weight_ext <- W_now[ind_range,] + weight_ext
    
    # Move one order forward
    W_cml <- W_cml_current + W_cml
    W_cml[W_cml != 0] <- 1
    
    # Control when graph is fully connected
    if(length(which(W %^% (i-1) > 1)) == n^2) break
  }
  
  return(weight_ext)
}

# Wrapper
f_mm_mat <- function(
    W, # Adjacency matrix
    m, # Number of memberships
    ord = 3, # Order of neighbours to include
    w_ord, # Weight of neighbours orders
    id_vec, # Whether an area is included in the membership (DEFAULT: ALL AREAS)
    excess_areas = F, # Areas to repeat if m > n
    red_areas, # Areas to not repeat if m < n
    rand_excess = F # Shuffle weights randomly for memberships
){
  library(expm)
  # browser()
  # Number of neighbours
  n <- nrow(W) ; if(n < 4) stop("Too few areas")
  
  ### Create weights for each new membership
  ### If the area identically included in the membership does not renormalise
  ### otherwise it will have to
  
  # General checks for weights
  ## First order (is area identically include in the membership?)
  if(exists("id_vec") == F) id_vec <- rep(1, n)
  ## Length of neighbour order weights vector
  if(length(w_ord) < (ord + 1)) stop("Not enough weights")
  ## Force 0 on weight for identical area weight if necessary
  if(sum(id_vec) == 0) w_ord[1] <- 0
  ## Check weights sum to 1
  if(sum(w_ord[1:(ord+1)]) != 1){
    warning(
      "Neighbour weights do not sum to one. Function will normalise and proceed"
    )
    w_ord <- abs(w_ord)/sum(w_ord)
  }
  
  # Identical area inclusion weight 
  weights_W1 <- diag(id_vec)*w_ord[1]
  
  # Normalise weights for membership without identical area
  ## Create register of weights for higher orders
  mat_ord_weight <- cbind(
    id_vec*w_ord[1], # Weights of identical order
    matrix(rep(w_ord[-1], n), nrow = n, ncol = length(w_ord[-1]), byrow = T)
  )
  ## Normalise register
  mat_ord_weight <- mat_ord_weight/rowSums(mat_ord_weight)
  
  # Compute necessary multiple order of neighbours matrices
  W_cml <- W
  weight_ext <- f_mm_high_ord(W, ord = ord, mat_ord_weight = mat_ord_weight,
                              ind_range = 1:n)
  
  # Add Identical neighbours and first order
  ## First order neigh
  weights_W2 <- W
  weights_W2[which(weights_W2 != 0)] <- runif(length(which(weights_W2 != 0)))
  for(j in 1:n){
    weights_W2[j, ] <- (weights_W2[j,]/sum(weights_W2[j,]))*mat_ord_weight[j,2]
  }
  ## Finalise
  weight_ext <- weight_ext + weights_W1 + weights_W2 ; rowSums(weight_ext)
  ## Check
  if(sum(rowSums(weight_ext)) != n){
    stop("Weights of weight_ext do not sum to 1")
  }
  
  # More MEMBERSHIPS than AREAS (m > n)
  if(m > n){
    # Check
    if(missing(excess_areas)){
      if(m - n > n){
        warning("m - n > n, so excess areas sampled with replacement")
        excess_areas <- sample(1:n, m - n, replace = T)
      } else{
        excess_areas <- sample(1:n, m - n)
      }
    } 
    
    # Create additional matrix
    ## Higher order
    excess_weight_high <- f_mm_high_ord(
      W, ord = ord, mat_ord_weight = mat_ord_weight, 
      ind_range = excess_areas)
    ## First order
    excess_weight_2 <- W
    excess_weight_2[which(excess_weight_2 != 0)] <-
      runif(length(which(excess_weight_2 != 0)))
    for(j in excess_areas){
      excess_weight_2[j, ] <- 
        (excess_weight_2[j,]/sum(excess_weight_2[j,]))*mat_ord_weight[j,2]
    }
    ## Subset to excess areas in second order
    excess_weight_2 <- excess_weight_2[excess_areas,]
    ## Identical areas
    excess_weight_1 <- matrix(0, nrow = length(excess_areas), ncol = n)
    excess_weight_1[cbind(1:length(excess_areas), excess_areas)] <- w_ord[1]
    ## Finalise
    excess_weight <- excess_weight_1 + excess_weight_2 + excess_weight_high
    
    # Final save
    weight_ext <- rbind(weight_ext, excess_weight)
    
    ## Check
    if(sum(rowSums(weight_ext)) != m){
      stop("Weights of weight_ext do not sum to 1")
    }
  }
  
  # More AREAS than MEMBERSHIPS (m < n)
  if(m < n){

    weight_ext <- weight_ext[sample(1:n, m), ]
  }
  
  # Return values
  output <- list(
    weight_ext = weight_ext,
    w_ord = w_ord,
    id_vec = id_vec,
    weights_W1 = weights_W1,
    weights_W2 = weights_W2,
    W_cml = W_cml
  )
  return(output)
}
#

#### CAR SIMULATION ####

# Univariate CAR
f_phi_car <- function(W, alpha = .5, tau = 10, B = 10){
  library(MASS)
  
  # Number of Neighbours per area
  D <- diag(rowSums(W))
  
  # Covariance matrix for CAR
  Q <- solve(tau*(D - alpha*W))
  
  # Generate phi's
  phi <- mvrnorm(n = B, mu = rep(0, nrow(W)), Sigma = Q)
  # colnames(phi) <- 1:nrow(W)
  
  return(phi)
}
#

#### SIMULATION BASED CALIBRATION ####

# Function for bias
f_bias <- function(x, true){
  c(mean(x - true), mean(abs(x - true)), sqrt(mean((x - true)^2)))
}
# Function for summary statistics
f_stat <- function(x){
  c(mean(x), sd(x), quantile(x, prob = c(.025, .05, .5, .95, .975)))
}

# Create array and matrix fits - works for m >= n
f_ar_mat <- function(
    fit, inv = F, niter, nchains, 
    pars_mean, phi_names, l_RR_names, rr_mm_names,
    H_inv = NULL
){
  # browser()
  
  n_adj <- length(l_RR_names) ; m <- length(rr_mm_names)
  
  # Parameter names
  all_pars <- c(pars_mean, phi_names, l_RR_names, rr_mm_names)
  
  # Create array and matrix
  ## Array
  ### INVERSE
  if(inv == T){
    # Create array
    array_inv <- as.array(fit)
    ## Names for array
    names_array_inv <- c(dimnames(array_inv)$parameters, l_RR_names)
    ## Create new array with appropriate parameters
    array_new <- array(
      NA, dim = c(niter/2, nchains, length(names_array_inv)),
      dimnames = list(
        iterations = NULL, 
        chains = sapply(1:nchains, function(x) paste("chain:", x, sep = "")),
        parameters = names_array_inv
      )
    )
    ## Copy already sampled values
    array_new[,, dimnames(array_inv)$parameters] <- 
      array_inv[,, dimnames(array_inv)$parameters]
    
    # Transform samples
    for(i in 1:nchains){
      ## Obtain samples for chain i
      chain <- apply(
        array_inv[,i, rr_mm_names[1:n_adj]], 1, function(x) H_inv %*% x[1:n_adj]
      ) %>% t()
      
      ## Save them in correct dimensions
      array_new[, i, l_RR_names] <- chain
      
      # Save
      array_fit <- array_new
    }
    
    # Matrix
    ## Compute inverted values
    l_RR_inv <- apply(
      as.matrix(fit)[, rr_mm_names], 1, function(x) H_inv %*% x[1:n]
    ) %>% t()
    ## Save them
    mat_fit <- matrix(NA, nrow = B, ncol = length(all_pars))
    colnames(mat_fit) <- all_pars
    mat_fit[, c(pars_mean, phi_names, rr_mm_names)] <- 
      as.matrix(fit)[, c(pars_mean, phi_names, rr_mm_names)]
    mat_fit[, l_RR_names] <- l_RR_inv
    
    ### POST
  } else{
    array_fit <- as.array(
      fit, pars = all_pars
    )
    mat_fit <- as.matrix(fit)
  }
  
  return(list(array = array_fit, mat = mat_fit))
}

# Compute SBC
f_sbc <- function(
    fit, index,
    # Names
    pars_mean, phi_names, l_RR_names, RR_names, rr_mm_names, 
    # Parameters
    list_pars, niter, nchains, thinning = 3,
    # Inverse
    H_inv = NULL, inv = F, unst = F,
    # Matrix and array
    array_fit, mat_fit
){
  # browser()
  # Libraries
  library(rstan)
  library(dplyr)
  
  # Parameters
  B <- nrow(as.matrix(fit)) ; n <- length(l_RR_names) ; m <- length(rr_mm_names)
  
  # Save results matrix
  all_pars <- c(pars_mean, phi_names, l_RR_names, RR_names, rr_mm_names)
  results <- matrix(NA, nrow = length(all_pars), ncol = 14)
  colnames(results) <- c("rank_stat", "Rhat", "ESS_bulk", "index",
                         "bias", "abs_bias", "RMSE",
                         "mean", "sd", "q025", "q05", "q5", "q95", "q975")
  rownames(results) <- all_pars
  
  # Fetch true parameters
  ## Mean
  true_pars <- list(
    alpha = list_pars$alpha, tau = list_pars$tau, gamma = list_pars$gamma, 
    `beta[1]` = list_pars$beta[, 1], `beta[2]` = list_pars$beta[, 2]
  )
  ## Areal 
  phi <- list_pars$phi ; l_RR <- d_grid$l_RR ; l_RRmm <- d_grid$l_RRmm
  ## Memberships
  
  # Assign index
  results[, "index"] <- rep(index, length(all_pars))
  ## Mean function parameters
  for(i in 1:length(pars_mean)){
    # Thinning
    ## ESS
    n_eff <- ess_bulk(array_fit[,,pars_mean[i]])
    ## Rank statistic on thinned
    current_rank <- sum(
      as.numeric(array_fit[,,pars_mean[i]])[seq(1, B, by = thinning)]
      < true_pars[[pars_mean[i]]][index]
    )
    
    # Save results
    ## Ranks and convergence
    results[pars_mean[i], c("rank_stat", "Rhat", "ESS_bulk")] <- c(
      current_rank, Rhat(array_fit[,,pars_mean[i]]), n_eff
    )
    ## Bias
    results[pars_mean[i], c("bias", "abs_bias", "RMSE")] <- f_bias(
      x = mat_fit[,pars_mean[i]], true = true_pars[[pars_mean[i]]][index]
    )
    ## Summary statistics
    results[pars_mean[i], c("mean", "sd", "q025", "q05", "q5", "q95", "q975")] <- 
      f_stat(x = mat_fit[,pars_mean[i]])
  }
  ## CAR
  for(i in 1:length(phi_names)){
    # Thinning
    ## ESS
    n_eff <- ess_bulk(array_fit[,,phi_names[i]])
    ## Rank statistic on thinned
    current_rank <- sum(
      as.numeric(array_fit[,,phi_names[i]])[seq(1, B, by = thinning)]
      < phi[index, i]
    )
    
    # Save results
    ## Ranks and convergence
    results[phi_names[i], c("rank_stat", "Rhat", "ESS_bulk")] <- c(
      current_rank, Rhat(array_fit[,,phi_names[i]]), n_eff
    )
    ## Bias
    results[phi_names[i], c("bias", "abs_bias", "RMSE")] <- f_bias(
      x = mat_fit[,phi_names[i]], true = phi[index, i]
    )
    ## Summary statistics
    results[phi_names[i], c("mean", "sd", "q025", "q05", "q5", "q95", "q975")] <- 
      f_stat(x = mat_fit[,phi_names[i]])
  }
  
  ## Areal RR
  for(i in 1:length(l_RR_names)){
    # Thinning
    ## ESS
    n_eff <- ess_bulk(array_fit[,,l_RR_names[i]])
    ## Rank statistic on thinned
    current_rank <- sum(
      as.numeric(array_fit[,,l_RR_names[i]])[seq(1, B, by = thinning)]
      < l_RR[index, i]
    )
    
    # Save results
    ## Ranks and convergence
    results[RR_names[i],  c("rank_stat", "Rhat", "ESS_bulk")] <- c(
      current_rank, Rhat(array_fit[,,l_RR_names[i]]), n_eff
    )
    ## Bias
    results[RR_names[i], c("bias", "abs_bias", "RMSE")] <- f_bias(
      x = exp(mat_fit[,l_RR_names[i]]),
      true = exp(l_RR[index, i])
    )
    ## Summary statistics
    results[RR_names[i], c("mean", "sd", "q025", "q05", "q5", "q95", "q975")] <- 
      f_stat(x = exp(mat_fit[,l_RR_names[i]]))
  }
  
  ## Membership RR
  for(i in 1:length(rr_mm_names)){
    # Thinning
    ## ESS
    n_eff <- ess_bulk(array_fit[,,rr_mm_names[i]])
    ## Rank statistic on thinned
    current_rank <- sum(
      as.numeric(array_fit[,,rr_mm_names[i]])[seq(1, B, by = thinning)]
      < l_RRmm[index, i]
    )
    
    # Save results
    ## Ranks and convergence
    results[rr_mm_names[i], c("rank_stat", "Rhat", "ESS_bulk")] <- c(
      current_rank, Rhat(array_fit[,,rr_mm_names[i]]), n_eff
    )
    ## Bias
    results[rr_mm_names[i], c("bias", "abs_bias", "RMSE")] <- f_bias(
      x = exp(mat_fit[,rr_mm_names[i]]),
      true = exp(l_RRmm[index, i])
    )
    ## Summary statistics
    results[rr_mm_names[i], c("mean", "sd", "q025", "q05", "q5", "q95", "q975")] <- 
      f_stat(x = exp(mat_fit[,rr_mm_names[i]]))
  }
  
  return(results)
}


# Functions
## Scaling functions
f_scale <- function(x) x
f_inv_scale <- function(x) x
## General conf bands
f_q_beta_low_p <- function(n, j, p){
  qbeta(p, shape1 = j*n, shape2 = n - j*n + 1)
}
f_q_beta_low_p <- Vectorize(f_q_beta_low_p, vectorize.args = "j")
f_q_beta_high_p <- function(n, j, p){
  qbeta(1-p, shape1 = j*n, shape2 = n - j*n + 1)
}
f_q_beta_high_p <- Vectorize(f_q_beta_high_p, vectorize.args = "j")
## SCALE - conf bands
f_scale_q_beta_low_p <- function(n, j, p){
  j <- f_inv_scale(j)
  f_scale(qbeta(p, shape1 = j*n, shape2 = n - j*n + 1))
}
f_scale_q_beta_low_p <- Vectorize(f_scale_q_beta_low_p, vectorize.args = "j")
f_scale_q_beta_high_p <- function(n, j, p){
  j <- f_inv_scale(j)
  f_scale(qbeta(1-p, shape1 = j*n, shape2 = n - j*n + 1))
}
f_scale_q_beta_high_p <- Vectorize(f_scale_q_beta_high_p, vectorize.args = "j")
## Coverage
f_coverage_count <- function(x, n, o_stat, p){
  cov_count <- 0
  for(j in 1:n){
    ifelse(
      x[j] > f_q_beta_high_p(n = n, j = o_stat[j], p = p) | 
        x[j] < f_q_beta_low_p(n = n, j = o_stat[j], p = p),
      cov_count <- cov_count + 1, cov_count <- cov_count + 0
    )
  }
  return(cov_count)
}

# PLOT SBC 
f_sbc_plot <- function(d_rank, v_pars, L_thin, gg_base, p_cov, vec_colours){
  # browser()
  # Obtain general parameters
  N <- nrow(d_rank)
  ord_stat <- 1:N/(N + 1)
  ## Number of parameters included in the plot
  D <- length(v_pars)
  
  # Confidence bands
  c_bands <- list(
    geom_function(fun = f_scale_q_beta_low_p, args = list(n = N+1, p = p_cov), 
                  col = "red", size = 1),
    geom_function(fun = f_scale_q_beta_high_p, args = list(n = N+1, p = p_cov), 
                  col = "red", size = 1)
  )
  
  # Loop for plot creation and coverage
  gg_sbc <- gg_base
  par_cov <- numeric()
  for(i in 1:D){
    current <- data.frame(
      ord = f_scale(1:N/(N + 1)), 
      r_vec = f_scale(sort(d_rank[,v_pars[i]])/L_thin))
    current_long <- current %>%
      pivot_longer(!ord, names_to = "names", values_to = "values") 
    ## Save plot
    gg_sbc <- gg_sbc +
      geom_point(data = current_long, aes(x = ord, y = values),
                 colour = "blue")
    ## Compute coverage
    par_cov[i] <- f_coverage_count(
      x = f_inv_scale(current$r_vec), n = N, o_stat = ord_stat, p = p_cov
    )
  }
  # Summarise coverage
  tot_coverage_perc <- sum(par_cov)/(D*N)*100
  return(gg_sbc + c_bands + annotate("text", x = .2, y = .9, size = 6, 
                                     label = paste("Coverage: ", 
                                                   round(100 - tot_coverage_perc, 3), "%"))
  )
}

# Base plot
base_plot <- ggplot() +
  # geom_function(fun = f_scale, col = "black", size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "none") + 
  xlab("") + ylab("") + xlim(0,1)
#

#### MIXED REPLICATES ####

# CAR-MM mixed replicates
f_mix_carmm <- function(alpha, tau, beta, gamma, covars, W, H, off, psi,
                        glm = F, icar = F, nb = T){
  # browser() 
  library(bayesSurv)
  # Diagonal matrix
  D_mat <- diag(rowSums(W))
  n <- nrow(W)
  
  # Bind vector for intercept
  covars <- cbind(rep(1, nrow(covars)), covars)
  
  # GLM CHECK
  if(glm == F){
    # Generate new phi CAR
    if(icar == F){
      phi_new <- t(apply(data.frame(alpha, tau), 1,
                         function(x) rMVNorm(n = 1, mean = 0, 
                                             Q = x[2]*(D_mat - x[1]*W),     
                                             param = "canonical")))
    } else{# Generate ICAR (martinez-beneito 2019 section 6.1)
      Q <- f_prec_car(W = W, alpha = 1, tau = 1)
      V <- eigen(Q)$vectors
      L <- diag(sqrt(1/(eigen(Q)$values)[1:(n-1)]))
      phi_new <- t(sapply(tau, function(x) 
        (1/sqrt(x))*V[,1:(n-1)] %*% L %*% rnorm(n-1)))
    }
  } else{
    phi_new <- matrix(0, nrow = length(gamma), ncol = nrow(W))
  }
  
  # Covariates effect
  covs_new <- t(apply(cbind(gamma, beta), 1, function(x) covars%*%x))
  
  # Generate RR
  l_RR_mix <-  covs_new + phi_new
  
  # Generate new outcomes
  if(nb == T){
    y_rep <- t(
      apply(
        cbind(psi, l_RR_mix), 1, 
        function(x) rnbinom(
          nrow(H), mu = off*exp(H %*% x[-1]), size = x[1]
        )
      )
    )
  } else{
    y_rep <- t(apply(
      l_RR_mix, 1, function(x) rpois(nrow(H), lambda = off*exp(H %*% x))
    ))
  }
  
  output <- list(
    y_rep = y_rep,
    RR_mix = exp(l_RR_mix),
    phi_mix = phi_new
  )
  return(output)
}
#


#### PPP - MARGINAL ####

# PPP discrete
f_ppp_disc <- function(y, yrep){
  if(nrow(yrep) != length(y)) stop("Dimensions of arguments do not match")
  apply(cbind(y, yrep), 1, 
        function(x) mean(x[-1] < x[1]) + .5*mean(x[-1] == x[1])
  )
}

# PPP continuous
f_ppp_cont <- function(y, yrep){
  if(nrow(yrep) != length(y)) stop("Dimensions of arguments do not match")
  apply(cbind(y, yrep), 1, function(x) mean(x[-1] < x[1]))
}

# Q-Q plot for marginal posterior checks
qq_marg_ppp <- function(p_vals){
  plot(sort(p_vals), (1:length(p_vals))/(length(p_vals) + 1))
  abline(a = 0, b = 1)
}
#

#### SCORING RULES ####

# Ranked Probability score
rps <- function(y_rep, y_obs){
  # Differences with observed values
  diff_real <- t(apply(y_rep, 1, function(x) abs(x - y_obs)))
  m_real <- apply(diff_real, 2, function(x) mean(abs(x)))
  
  # Floor of half of simulations
  n_2 <- floor(nrow(y_rep)/2)
  
  # Difference among generated repetitions
  half_diff_rep <- matrix(NA, ncol = ncol(y_rep), nrow = n_2)
  for(j in 1:ncol(y_rep)){
    for(i in 1:n_2){
      half_diff_rep[i, j] <- abs(y_rep[i,j] -  y_rep[i+n_2,j])
    }
  }
  m_rep <- t(apply(half_diff_rep, 2, function(x) sum(x)/nrow(y_rep)))
  
  mean(m_real - m_rep)
}

# Dawid-Sebastiani score
dss <- function(y_rep, y_obs){
  m_rep <- t(apply(y_rep, 2, mean))
  sd_rep <- t(apply(y_rep, 2, sd))
  
  mean(((y_obs - m_rep)/sd_rep)^2 + 2*log(sd_rep))
}
#