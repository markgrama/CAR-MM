#### MM MATRIX SIMULATION ####

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
  
  # Create matrix of weights
  # weight_ext <- matrix(0, ncol = n, nrow = n)
  
  ### Create weights for each new membership
  ### If the area identically included in the membership does not renormalise
  ### otherwise it will have to
  
  # General checks for weights
  ## First order (is area identically include in the membership?)
  # if(exists("id_vec") == F) id_vec <- rep(1, n)
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
    
    #### DEPRECATED
    # Shuffle excess memberships weights if required
    # if(rand_excess == T){
    #   for(i in 1:(m-n)){
    #     excess_weight[i] <- weight_ext[excess_areas[i], sample(1:n, n)]
    #   }
    # } else{
    #   excess_weight <- weight_ext[excess_areas, ]
    # }
    
    # Final save
    weight_ext <- rbind(weight_ext, excess_weight)
    
    ## Check
    if(sum(rowSums(weight_ext)) != m){
      stop("Weights of weight_ext do not sum to 1")
    }
  }
  
  # More AREAS than MEMBERSHIPS (m < n)
  if(m < n){
    # if(rand_excess == T){
    #   
    #   # Create additional matrix
    #   weight_ext_f <- matrix(0, ncol = n, nrow = m - n)
    #   
    #   for(i in 1:m){
    #     weight_ext_f[i] <- weight_ext[red_areas[i], sample(1:n, n)]
    #   }
    #   # Save final matrix
    #   ## Check
    #   if(sum(rowSums(weight_ext)) != n){
    #     stop("Weights of weight_ext do not sum to 1")
    #   }
    #   weight_ext <- weight_ext_f
    # } else{
    #   # Save final matrix
    #   weight_ext <- weight_ext[red_areas, ]
    #   ## Check
    #   if(sum(rowSums(weight_ext)) != n){
    #     stop("Weights of weight_ext do not sum to 1")
    #   }
    # }
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

# Reverse MM matrix
f_reverse_mm_mat <- function(H, counts){
  # Control
  if(sum(round(rowSums(H), 0)) != nrow(H)) stop("Not an MM matrix")
  
  # Retrieve orignal matrix of counts
  og_mat <- H*counts
  
  # Output the row-normalised matrix for areas to GP
  t(og_mat)/rowSums(t(og_mat))
}

# Allocate remainder
f_allocate_remainder <- function(
  tot, # Quantity to allocate almost equally
  size # Number of dimensions to allocate to 
){
  # Compute remainder
  remainder <- tot - floor(tot / ceiling(tot/size))*ceiling(tot/size)
  ## Create vector
  vec_rem <- c(rep(ceiling(tot/size), floor(tot/ceiling(tot/size))), remainder)
  
  # Correct size
  ifelse(length(vec_rem) < size,
         vec <- c(vec_rem, rep(0,  size - length(vec_rem))),
         vec <- vec_rem)
  
  return(vec)
}

# Allocate seats
f_allocate_seats <- function(n, shares){
  # Price of a single seat
  price <- 1/n
  
  # Initial allocation
  seats <- floor(shares/price)
  ## Compute remainder
  remainder <- n - sum(seats)
  
  # Allocate remainder
  vec_rem <- f_allocate_remainder(tot = remainder, size = length(shares))
  # ifelse(length(vec_rem) < length(shares),
  # vec <- c(vec_rem, rep(0,  length(shares) - length(vec_rem))),
  # vec <- vec_rem)
  seats <- seats + vec_rem
  
  return(seats)
}

# Embed allocation of remainder
f_embed_rem <- function(vector, target_sum){
  # Check
  if(sum(vector) < target_sum){
    stop("Target sum under vector, use ceiling function instead")
  }
  
  as.integer(
    vector - f_allocate_remainder(tot = sum(vector) - target_sum,
                                  size = length(vector))
  )
}

# Full population simulation
f_full_pop_sim <- function(area_pop, rate_pop, out_rate,
                           H_inv_init, H_init,
                           W){
  ## Number of memberships and areas
  m <- nrow(H_init) ; n <- ncol(H_init)
  # Check for membership matrix
  if(dim(H_init)[2] != nrow(W)){
    stop("MM (inv) matrix does not match number of areas")
  }
  
  # Assign area population to memberships
  pop_memb <- ceiling(t(H_inv_init) %*% area_pop)
  
  # Assign area population to memberships
  pop_memb <- ceiling(t(H_inv_init) %*% area_pop)
  
  # Adjust populations
  ## Difference
  diff_pop <- sum(pop_memb) - sum(area_pop)
  ## Remove extra individuals
  pop_memb <- pop_memb - 
    f_allocate_remainder(tot = diff_pop, size = m)
  ## Check
  if(sum(pop_memb) != sum(area_pop)){
    stop("Membership and area populations do not match")
  }
  
  # Create structure population
  ## Area
  pop_str_area <- t(sapply(
    area_pop, function(x) f_allocate_seats(n = x, shares = rate_pop)
  ))
  ## Membership
  pop_str_memb <- t(sapply(
    pop_memb, function(x) f_allocate_seats(n = x, shares = rate_pop)
  ))
  
  # Compute offsets
  off_areas <- apply(pop_str_area, 1, function(x) x%*%out_rate)
  off_memb <- apply(pop_str_memb, 1, function(x) x%*%out_rate)
  
  # Recompute H accordingly
  H_abs <- matrix(0, nrow = m, ncol = n)
  for(i in 1:nrow(H_abs)){
    H_abs[i,] <- f_allocate_seats(
      n = pop_memb[i], shares = H_init[i,]
    )
  }
  ## Final matrices
  H <- H_abs/rowSums(H_abs)
  H_inv <- t(apply(H_abs, 2, function(x) x/sum(x)))
  
  # Bind output
  output <- list(
    pop_str_area = pop_str_area,
    pop_str_memb = pop_str_memb,
    off_areas = off_areas,
    off_memb = off_memb,
    H = H,
    H_inv = H_inv,
    H_abs = H_abs
  )
  return(output)
}
#

#### DAGAR ####

# Q matrix for DAGAR
q_DAGAR <- function(W, rho, bivariate = F){
  # Control for disconnected graph
  if(sum(rowSums(W) == 0) >= 1) print("Error: disconnected graph")
  
  # Number of areas
  n <- nrow(W)
  
  # Retrieve univariate rho
  rho1 <- rho[1]
  
  # Creation of B_adj which is going to be the adj matrix of the graph
  B_adj <- W ; B_adj[upper.tri(B_adj)] <- 0
  
  B1 <- B_adj
  for(i in 1:nrow(B1)){
    B1[i, which(B1[i,] == 1)] <- rho1/(1 + (sum(B1[i,]) - 1)*rho1^2)
  }
  
  # Compute precision matrix for w
  ## Obtain L
  L1 <- diag(rep(1, n)) - B1
  ## Compute precisions
  tau1 <- (1 + (apply(B_adj, 1, sum) - 1)*rho1^2)/(1 - rho1^2)
  
  if(bivariate == T){
    # Retrieve 2nd rho
    rho2 <- rho[2]
    
    B2 <- B_adj
    for(i in 1:nrow(B2)){
      B2[i, which(B2[i,] == 1)] <- rho2/(1 + (sum(B2[i,]) - 1)*rho2^2)
    }
    
    ## Obtain L
    L2 <- diag(rep(1, n)) - B2
    ## Compute precisions
    tau2 <- (1 + (apply(B_adj, 1, sum) - 1)*rho2^2)/(1 - rho2^2)
    
    # Save output
    output <- list(Q1 = t(L1)%*%diag(tau1)%*%L1,
                   Q2 = t(L2)%*%diag(tau2)%*%L2)
    return(output)
  } else{
    Q <- t(L1)%*%diag(tau1)%*%L1 ; return(Q) 
  }
}

# Order free DAGAR covariance (from Datta script)
q_of_DAGAR <- function(W, rho){
  # Parameter
  u <- rho^2
  
  # Number of neighbours for each area
  Minc <- W
  ni <- rowSums(Minc)
  n <- length(ni)
  maxn <- max(ni)
  
  # Equivalent of -f_fRho-
  sumfuncvec <- numeric()
  for(i in 1:maxn)  sumfuncvec[i]=i/(1-u+i*u)
  
  # Cumulative sum for diagonal terms
  cumsumfuncvec <- numeric()
  cumsumfuncvec[1] = sumfuncvec[1]
  for(i in 2:maxn) cumsumfuncvec[i] = cumsumfuncvec[i-1]+sumfuncvec[i]
  
  # Diagonal matrix    
  Qd <- matrix(0, n, n)
  
  # Neighbours list
  neighbors <- lapply(1:n, function(x) which(Minc[x,] == 1))
  
  # Base quantities
  neimat <- matrix(0, n, maxn) # considers the largest no off neighbours
  for(i in 1:n) neimat[i,1:ni[i]] <- neighbors[[i]]
  ##
  intersectmat <- matrix(0,n^2,maxn)
  nijvec <- rep(0,n^2)
  for(i in 1:n) for(j in 1:n)
  { neighij=intersect(neighbors[[i]],neighbors[[j]])
  nij=length(neighij)
  nijvec[n*(i-1)+j]=nij
  if(nij >0) intersectmat[n*(i-1)+j,1:nij]=neighij
  }
  
  for(i in 1:n){
    s <- 0
    if(ni[i]>0) for(ck in 1:ni[i]){
      # Neimat is the index lookup matrix for areas neighbours
      k = neimat[i,ck]
      # Summation term in diagonal matrix with rho^2
      s = s + u*cumsumfuncvec[ni[k]]/(ni[k]*(ni[k]+1))
      # Common term for off-diagonal elements
      Qd[i,k] = Qd[i,k] - rho
    }
    
    # Almost finished diagonal term
    Qd[i,i] = 1 - u + u*ni[i]/2 + s
    
    if(i < n){
      for(j in (i+1):n){
        t=0
        jc=0
        counter=(i-1)*n+j
        nij=nijvec[counter]
        if(nij>0){
          jointn=intersectmat[counter,1:nij]
          for(ck in 1:length(jointn))
          {
            k=jointn[ck]
            t=t+1/(2*(ni[k]+1))+(ni[k]-cumsumfuncvec[ni[k]])/(ni[k]*(ni[k]+1)*(ni[k]-1))
          }
        }
        Qd[i,j]=Qd[i,j]+t
        Qd[j,i]=Qd[j,i]+t
      }
    }
  }    
  
  # Final common term for diagonal term
  Qd <- Qd/(1-rho^2)
  return(Qd)
}

# Generation DAGAR random effects
f_phi_dagar <- function(W, rho, tau, eta, B = 1E3, bivariate = F, OF = F){
  library(MASS)
  
  # Number of areas
  n <- nrow(W)
  
  # Bivariate
  if(bivariate == T){
    # Choose which DAGAR prior
    if(OF == T){
      # ORDER-FREE
      Q_biv <- list(Q1 = q_of_DAGAR(W, rho[1]),
                    Q2 = q_of_DAGAR(W, rho[2]))
    } else{
      # ORDERED
      Q_biv <- q_DAGAR(W, rho, bivariate = T)
    }
    
    # Etas
    eta_0 <- eta[1] ; eta_1 <- eta[2]
    
    # Phi
    phi <- matrix(NA, ncol = 2*n, nrow = B)
    ## Phi 2
    phi[, (n + 1):(2*n)] <- mvrnorm(n = 1E3, mu = rep(0, n), 
                                    Sigma = solve(Q_biv$Q2*tau[2]))
    phi2 <- phi[, (n + 1):(2*n)]
    ## Phi 1
    phi[, 1:n] <- mvrnorm(n = 1E3, mu = rep(0, n), 
                          Sigma = solve(Q_biv$Q1*tau[1]))
    phi1 <- phi[, 1:n] + 
      t(apply(phi2, 1, function(x) (eta_0*diag(rep(1, n)) + eta_1*W_eW)%*%x))
    
    return(list(phi1 = phi1, phi2 = phi2, Q = Q_biv))
  } else{
    # Q matrix
    if(OF == T){
      Q <- q_of_DAGAR(W, rho)
    } else{
      Q <- q_DAGAR(W, rho)
    }
    
    # Phi
    phi <- mvrnorm(n = B, mu = rep(0, n), Sigma = solve(Q*tau[1]))
    
    return(list(phi = phi, Q = Q))
  }
}

# Outcome
f_dagar_out <- function(W, M_W, rho, eta, off1, off2, 
                        B = 1E3, index = c(1,1), overd = 10){
  # Number of observations
  n <- nrow(W) # Areas
  m <- nrow(M_W) # Memberships
  
  # Phi
  dagar <- f_phi_dagar(W, rho = rho, eta = eta, B = B, bivariate = T)
  phi1 <- dagar$phi1
  phi2 <- dagar$phi2
  
  # Generate outcomes
  y1 <- rnbinom(n, mu = off1*exp(phi1[index[1],]), size = overd)
  y2 <- rnbinom(m, mu = off2*exp(M_W%*%phi2[index[2],]), size = overd)
  
  return(list(
    Q = dagar$Q_biv,
    phi1 = phi1, phi2 = phi2,
    y1 = y1, y2 = y2,
    off1 = off1, off2 = off2
  ))
}

# Average neighbour correlation
f_avgNegCor <- function(W, phi){
  l_cor <- list()
  for(i in 2:(nrow(W)-1)){
    if(length(which(W[i,] == 1)) == 1) next
    l_cor[[i]] <- apply(phi[, which(W[i,] == 1)], 2, 
                        function(x) mean(cor(x, phi[,i])))
  }
  mean(unlist(l_cor), na.rm = T)
}

# Adjacency objects DAGAR-OF
f_adj_dagar_of <- function(W){
  n <- nrow(W)
  n_i <- rowSums(W)
  maxn_i <- max(n_i)
  # NEIMAT
  ## Neighbours list
  neighbors <- lapply(1:n, function(x) which(W[x,] == 1))
  ## Base quantities
  neimat <- matrix(0, n, maxn_i) # considers the largest no off neighbours
  for(i in 1:n) neimat[i, 1:n_i[i]] <- neighbors[[i]]
  
  # N_IJVEC
  intersectmat <- matrix(0, n^2, maxn_i)
  n_ijvec <- rep(0, n^2)
  for(i in 1:n) for(j in 1:n){
    neigh_ij <- intersect(neighbors[[i]], neighbors[[j]])
    n_ij <- length(neigh_ij)
    n_ijvec[n*(i-1) + j] <- n_ij
    if(n_ij >0) intersectmat[n*(i-1) + j, 1:n_ij] <- neigh_ij
  }
  
  output <- list(n = n, n2 = n^2, n_i = n_i,
                 maxn_i = maxn_i, neimat = neimat,
                 n_ijvec = n_ijvec, intersectmat = intersectmat)
  return(output)
}

# Stan OF-DAGAR Data simulation
f_stan_data <- function(
  # Mandatory
  W, log_offset1, rho, tau, gamma_i,
  # Optional
  log_offset2, X_cov, eta, beta1, beta2, overd,
  # Controls
  cov = F, biv = F, type = "POI", OF = T){
  
  # Common elements
  n <- nrow(W)
  n_i <- rowSums(W)
  maxn_i <- max(n_i)
  # NEIMAT
  ## Neighbours list
  neighbors <- lapply(1:n, function(x) which(W[x,] == 1))
  ## Base quantities
  neimat <- matrix(0, n, maxn_i) # considers the largest no off neighbours
  for(i in 1:n) neimat[i, 1:n_i[i]] <- neighbors[[i]]
  
  # N_IJVEC
  intersectmat <- matrix(0, n^2, maxn_i)
  n_ijvec <- rep(0, n^2)
  for(i in 1:n) for(j in 1:n){
    neigh_ij <- intersect(neighbors[[i]], neighbors[[j]])
    n_ij <- length(neigh_ij)
    n_ijvec[n*(i-1) + j] <- n_ij
    if(n_ij >0) intersectmat[n*(i-1) + j, 1:n_ij] <- neigh_ij
  }
  
  # Save common elements
  adj_base_list <- list(
    W = W,
    neimat = neimat,
    n2 = n^2,
    n = n,
    n_i = n_i, 
    maxn_i = maxn_i,
    n_ijvec = n_ijvec,
    intersectmat = intersectmat
  )
  
  # Generate relative risk (RR)
  ## UNIVARIATE
  if(biv == F){
    # Random effects
    if(OF == T){
      phi_dagar <- f_phi_dagar(W, rho = rho, tau = tau, B = 1,
                               bivariate = F, OF = T)$phi
    }
    if(OF == F){
      phi_dagar <- f_phi_dagar(W, rho = rho, tau = tau, B = 1,
                               bivariate = F, OF = F)$phi
    }
    
    # Covariate effect
    if(cov == T){
      X_i <- X_cov%*%beta1 
      k <- ncol(X_cov)
    } else{
      X_i <- rep(0, n)
      X_cov <- rep(0, n)
      k <- 0
    }
    
    # RR
    l_RR <- phi_dagar + gamma_i + X_i
    
    # Outcome simulation
    if(type == "POI"){ 
      y1 <- rpois(n, exp(log_offset1 + l_RR))
      ## Save
      par_l <- list(beta1 = beta1, gamma_i = gamma_i,
                    phi = phi_dagar, k = k, tau = tau, rho = rho)
      out_l <- list(y1 = y1, X_cov = X_cov, log_offset1 = log_offset1)
    }
    if(type == "NB"){
      y1 <- rnbinom(n, mu = exp(log_offset1 + l_RR), size = overd[1])
      ## Save
      par_l <- list(beta1 = beta1, gamma_i = gamma_i, overd = overd,
                    phi = phi_dagar, k = k, tau = tau, rho = rho)
      out_l <- list(y1 = y1, X_cov = X_cov, log_offset1 = log_offset1)
    }
    
    # Save output
    outcome <- append(out_l, par_l)
  } else{ ## BIVARIATE
    # Random effects
    if(OF == T){
      phi_dagarT <- f_phi_dagar(W, rho = rho, tau = tau, eta = eta,
                                B = 1, bivariate = T, OF = T)
      phi_dagar1 <- phi_dagarT$phi1 ; phi_dagar2 <- phi_dagarT$phi2
    }
    if(OF == F){
      phi_dagarT <- f_phi_dagar(W, rho = rho, tau = tau, eta = eta,
                                B = 1, bivariate = T, OF = F)
      phi_dagar1 <- phi_dagarT$phi1 ; phi_dagar2 <- phi_dagarT$phi2
    }
    
    # Covariate effect
    if(cov == T){
      X_i1 <- X_cov%*%beta1  ; X_i2 <- X_cov%*%beta2 
      k <- ncol(X_cov)
    } else{
      X_i1 <- rep(0, n) ; X_i2 <- rep(0, n)
      X_cov <- rep(0, n)
      k <- 0
    }
    
    # RR
    l_RR1 <- phi_dagar1 + gamma_i[1] + X_i1
    l_RR2 <- phi_dagar2 + gamma_i[2] + X_i2
    
    # Outcome simulation
    if(type == "POI"){
      y1 <- rpois(n, exp(log_offset1 + l_RR1))
      y2 <- rpois(n, exp(log_offset2 + l_RR2))
      ## Save
      par_l <- list(beta1 = beta1, beta2 = beta2, gamma_i = gamma_i,
                    phi1 = phi_dagar1, phi2 = phi_dagar2, k = k,
                    tau = tau, rho = rho)
      out_l <- list(y1 = y1, y2 = y2, X_cov = X_cov,
                    log_offset1 = log_offset1, log_offset2 = log_offset2)
    }
    if(type == "NB"){
      y1 <- rnbinom(n, mu = exp(log_offset1 + l_RR1), size = overd[1])
      y2 <- rnbinom(n, mu = exp(log_offset2 + l_RR2), size = overd[2])
      ## Save
      par_l <- list(beta1 = beta1, beta2 = beta2, gamma_i = gamma_i, 
                    overd = overd, phi1 = phi_dagar1, phi2 = phi_dagar2,
                    k = k, tau = tau, rho = rho)
      out_l <- list(y1 = y1, y2 = y2, X_cov = X_cov,
                    log_offset1 = log_offset1, log_offset2 = log_offset2)
    }
    
    # Save output
    outcome <- append(out_l, par_l)
  }
  
  output <- append(adj_base_list, outcome)
  return(output)
}
#

#### DAGARMM ####

# Stan OF-DAGAR Data simulation
f_stan_data_mm <- function(
  # Mandatory
  W, log_offset1, rho, tau, gamma_i, M_W,
  # Optional
  log_offset2, X_cov, eta, beta1, beta2, overd,
  # Controls
  cov = F, biv = F, type = "POI", OF = T){
  
  # Common elements
  n <- nrow(W)
  n_i <- rowSums(W)
  maxn_i <- max(n_i)
  # NEIMAT
  ## Neighbours list
  neighbors <- lapply(1:n, function(x) which(W[x,] == 1))
  ## Base quantities
  neimat <- matrix(0, n, maxn_i) # considers the largest no off neighbours
  for(i in 1:n) neimat[i, 1:n_i[i]] <- neighbors[[i]]
  
  # N_IJVEC
  intersectmat <- matrix(0, n^2, maxn_i)
  n_ijvec <- rep(0, n^2)
  for(i in 1:n) for(j in 1:n){
    neigh_ij <- intersect(neighbors[[i]], neighbors[[j]])
    n_ij <- length(neigh_ij)
    n_ijvec[n*(i-1) + j] <- n_ij
    if(n_ij >0) intersectmat[n*(i-1) + j, 1:n_ij] <- neigh_ij
  }
  
  # Save common elements
  adj_base_list <- list(
    W = W,
    neimat = neimat,
    n2 = n^2,
    n = n,
    n_i = n_i, 
    maxn_i = maxn_i,
    n_ijvec = n_ijvec,
    intersectmat = intersectmat
  )
  
  # Generate relative risk (RR)
  ## UNIVARIATE
  if(biv == F){
    # Random effects
    if(OF == T){
      phi_dagar <- f_phi_dagar(W, rho = rho, tau = tau, B = 1,
                               bivariate = F, OF = T)$phi
    }
    if(OF == F){
      phi_dagar <- f_phi_dagar(W, rho = rho, tau = tau, B = 1,
                               bivariate = F, OF = F)$phi
    }
    
    # Covariate effect
    if(cov == T){
      X_i <- X_cov%*%beta1 
      k <- ncol(X_cov)
    } else{
      X_i <- rep(0, n)
      X_cov <- rep(0, n)
      k <- 0
    }
    
    # RR
    l_RR <- M_W%*%(phi_dagar + gamma_i + X_i)
    
    # Outcome simulation
    if(type == "POI"){ 
      y1 <- rpois(nrow(M_W), exp(log_offset1 + l_RR))
      ## Save
      par_l <- list(beta1 = beta1, gamma_i = gamma_i,
                    phi = phi_dagar, k = k, m = nrow(M_W))
      out_l <- list(y1 = y1, X_cov = X_cov, log_offset1 = log_offset1,
                    M_W = M_W)
    }
    if(type == "NB"){
      y1 <- rnbinom(nrow(M_W), mu = exp(log_offset1 + l_RR), size = overd[1])
      ## Save
      par_l <- list(beta1 = beta1, gamma_i = gamma_i, overd = overd,
                    phi = phi_dagar, k = k, m = nrow(M_W))
      out_l <- list(y1 = y1, X_cov = X_cov, log_offset1 = log_offset1,
                    M_W = M_W)
    }
    
    # Save output
    outcome <- append(out_l, par_l)
  } else{ ## BIVARIATE
    # Random effects
    if(OF == T){
      phi_dagarT <- f_phi_dagar(W, rho = rho, tau = tau, eta = eta,
                                B = 1, bivariate = T, OF = T)
      phi_dagar1 <- phi_dagarT$phi1 ; phi_dagar2 <- phi_dagarT$phi2
    }
    if(OF == F){
      phi_dagarT <- f_phi_dagar(W, rho = rho, tau = tau, eta = eta,
                                B = 1, bivariate = T, OF = F)
      phi_dagar1 <- phi_dagarT$phi1 ; phi_dagar2 <- phi_dagarT$phi2
    }
    
    # Covariate effect
    if(cov == T){
      X_i1 <- X_cov%*%beta1  ; X_i2 <- X_cov%*%beta2 
      k <- ncol(X_cov)
    } else{
      X_i1 <- rep(0, n) ; X_i2 <- rep(0, n)
      X_cov <- rep(0, n)
      k <- 0
    }
    
    # RR
    l_RR1 <- phi_dagar1 + gamma_i[1] + X_i1
    l_RR2 <- M_W%*%(phi_dagar2 + gamma_i[2] + X_i2)
    
    # Outcome simulation
    if(type == "POI"){
      y1 <- rpois(n, exp(log_offset1 + l_RR1))
      y2 <- rpois(nrow(M_W), exp(log_offset2 + l_RR2))
      ## Save
      par_l <- list(beta1 = beta1, beta2 = beta2, gamma_i = gamma_i,
                    phi1 = phi_dagar1, phi2 = phi_dagar2, k = k, m = nrow(M_W))
      out_l <- list(y1 = y1, y2 = y2, X_cov = X_cov, M_W = M_W,
                    log_offset1 = log_offset1, log_offset2 = log_offset2)
    }
    if(type == "NB"){
      y1 <- rnbinom(n, mu = exp(log_offset1 + l_RR1), size = overd[1])
      y2 <- rnbinom(n, mu = exp(log_offset2 + l_RR2), size = overd[2])
      ## Save
      par_l <- list(beta1 = beta1, beta2 = beta2, gamma_i = gamma_i, overd = overd,
                    phi1 = phi_dagar1, phi2 = phi_dagar2, k = k, m = nrow(M_W))
      out_l <- list(y1 = y1, y2 = y2, X_cov = X_cov, M_W = M_W,
                    log_offset1 = log_offset1, log_offset2 = log_offset2)
    }
    
    # Save output
    outcome <- append(out_l, par_l)
  }
  
  output <- append(adj_base_list, outcome)
  return(output)
}
#

#### ADJACENCY ####

# Rotation matrix (counterclockwise if premultiply)
Rmat <- function(angle){
  matrix(c(cos(angle), -sin(angle),
           sin(angle), cos(angle)), ncol = 2, byrow = T)
}

# Grid adjacency
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

# Function to create adj matrix from T and F
adj_inferred <- function(W_check){
  W <- matrix(0, nrow = nrow(W_check), ncol = nrow(W_check))
  
  for(i in 1:nrow(W_check)){
    for(j in 1:nrow(W_check)){
      W[i, j] <- ifelse(W_check[i, j], 1, 0)
    }
  }
  diag(W) <- rep(0, nrow(W_check))
  
  return(W)
}
#

#### CAR ####

# Jin et al. 2005 univariate CAR
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

# ICAR
f_phi_intrinsic <- function(W, tau = 10, B = 1){
  n <- nrow(W)
  
  Q <- f_prec_car(W = W, alpha = 1, tau = 1)
  
  V <- eigen(Q)$vectors
  L <- diag(1/eigen(Q)$values[1:(n-1)])
  
  mat_phi <- matrix(NA, nrow = B, ncol = n)
  for(i in 1:B) mat_phi[i,] <- 
    sqrt(1/tau)*V[,1:(n-1)]%*%sqrt(L)%*%rnorm(n-1)
  
  return(mat_phi)
}

# Generate entire CAR model
f_stan_data_car_mm <- function(
  # Mandatory
  W, log_offset1, alpha, tau, gamma_i, M_W,
  # Optional
  log_offset2, X_cov, eta, beta1, beta2, overd,
  # Controls
  cov = F, biv = F, type = "POI"){
  
  # browser()
  n <- nrow(W)
  
  # Generate relative risk (RR)
  ## UNIVARIATE
  if(biv == F){
    # Random effects
    phi_car <- f_phi_car(W, alpha = alpha, tau = tau, B = 1)
    
    # Covariate effect
    if(cov == T){
      X_i <- X_cov%*%beta1 
      k <- ncol(X_cov)
    } else{
      X_i <- rep(0, n)
      X_cov <- rep(0, n)
      k <- 0
    }
    
    # RR
    l_RR <- M_W%*%(phi_car + gamma_i + X_i)
    
    # Outcome simulation
    if(type == "POI"){ 
      y1 <- rpois(nrow(M_W), exp(log_offset1 + l_RR))
      ## Save
      par_l <- list(beta1 = beta1, gamma_i = gamma_i,
                    phi = phi_car, k = k, m = nrow(M_W),
                    n = nrow(W), W_n = sum(W)/2, W = W)
      out_l <- list(y1 = y1, X_cov = X_cov, log_offset1 = log_offset1,
                    M_W = M_W)
    }
    if(type == "NB"){
      y1 <- rnbinom(nrow(M_W), mu = exp(log_offset1 + l_RR), size = overd[1])
      ## Save
      par_l <- list(beta1 = beta1, gamma_i = gamma_i, overd = overd,
                    phi = phi_car, k = k, m = nrow(M_W),
                    n = nrow(W), W_n = sum(W)/2, W = W)
      out_l <- list(y1 = y1, X_cov = X_cov, log_offset1 = log_offset1,
                    M_W = M_W)
    }
    
    # Save output
    outcome <- append(out_l, par_l)
  } else{ ## BIVARIATE
    # Random effects
    phi_carT <- f_phi_car(W, alpha = alpha, tau = tau, eta = eta,
                          B = 1)
    phi_car1 <- phi_carT ; phi_car2 <- phi_carT
    
    # Covariate effect
    if(cov == T){
      X_i1 <- X_cov%*%beta1  ; X_i2 <- X_cov%*%beta2 
      k <- ncol(X_cov)
    } else{
      X_i1 <- rep(0, n) ; X_i2 <- rep(0, n)
      X_cov <- rep(0, n)
      k <- 0
    }
    
    # RR
    l_RR1 <- phi_car1 + gamma_i[1] + X_i1
    l_RR2 <- M_W%*%(phi_car2 + gamma_i[2] + X_i2)
    
    # Outcome simulation
    if(type == "POI"){
      y1 <- rpois(n, exp(log_offset1 + l_RR1))
      y2 <- rpois(nrow(M_W), exp(log_offset2 + l_RR2))
      ## Save
      par_l <- list(beta1 = beta1, beta2 = beta2, gamma_i = gamma_i,
                    phi1 = phi_car1, phi2 = phi_car2, k = k, m = nrow(M_W),
                    n = nrow(W), W_n = sum(W)/2, W = W)
      out_l <- list(y1 = y1, y2 = y2, X_cov = X_cov, M_W = M_W,
                    log_offset1 = log_offset1, log_offset2 = log_offset2)
    }
    if(type == "NB"){
      y1 <- rnbinom(n, mu = exp(log_offset1 + l_RR1), size = overd[1])
      y2 <- rnbinom(n, mu = exp(log_offset2 + l_RR2), size = overd[2])
      ## Save
      par_l <- list(beta1 = beta1, beta2 = beta2, gamma_i = gamma_i, overd = overd,
                    phi1 = phi_car1, phi2 = phi_car2, k = k, m = nrow(M_W),
                    n = nrow(W), W_n = sum(W)/2, W = W)
      out_l <- list(y1 = y1, y2 = y2, X_cov = X_cov, M_W = M_W,
                    log_offset1 = log_offset1, log_offset2 = log_offset2)
    }
    
    # Save output
    outcome <- append(out_l, par_l)
  }
  
  output <- outcome
  return(output)
}

# Precision matrix for CAR
f_prec_car <- function(W, alpha = .5, tau = 10){
  # Number of Neighbours per area
  D <- diag(rowSums(W))
  
  # Covariance matrix for CAR
  tau*(D - alpha*W)
}

# Leroux
f_prec_leroux <- function(W, tau, lambda){
  D <- diag(rowSums(W))
  
  tau*(lambda*(D - W) + diag(rep((1 - lambda), nrow(W))))
}
#

#### GMCAR ####

# GMCAR-MM SIM
sim_data_Gmcar <- function(
  pars, adj, offs, m_memb = F, mcar = F, M_W = NA
){
  library(MASS) # Simulate multivariate normal
  
  # browser()
  # Parameter names
  ## alpha1 - alpha2 - eta0 - eta1 - tau1 - tau2 - 
  ## cov - beta1[2] - beta2[2] - int[2] - overd1 - 
  ## overd 2
  
  # Revert to MCAR from GMCAR
  if(mcar == T){pars$eta1 <- 0 ; pars$alpha2 <- pars$alpha1}
  
  # Covariance matrix (see Jin et. al 2005 paper)
  mat_reg <- pars$eta0*diag(1, adj$n) + pars$eta1*adj$W
  mat_phi1 <- solve(pars$tau1*(adj$D-pars$alpha1*adj$W))
  mat_phi2 <- solve(pars$tau2*(adj$D-pars$alpha2*adj$W))
  sigma_11 <- mat_phi1 + mat_reg%*%mat_phi2%*%mat_reg
  sigma_12 <- mat_reg%*%mat_phi2
  sigma_21 <- mat_phi2%*%mat_reg
  sigma <- rbind(cbind(sigma_11, sigma_12),
                 cbind(sigma_21, mat_phi2))
  
  # Spatial effects
  phi_tot <- mvrnorm(n = 1,
                     mu = c(rep(0, adj$n),
                            rep(0, adj$n)), 
                     Sigma = sigma)
  phi_1 <- phi_tot[1:adj$n]
  phi_2 <- phi_tot[(adj$n + 1):(2*adj$n)]
  
  # Multiple membership
  if(m_memb == T & is.na(sum(M_W)) == T){
    m <- pars$m_memb # number of misaligned areas
    M_W1 <- matrix(runif(adj$n*m), ncol = adj$n, nrow = m) # weights
    M_W1 <- M_W1*matrix(sample(0:1, adj$n*m, replace = T),
                        ncol = adj$n)
    M_W <- M_W1/rowSums(M_W1) # normalized weights
    append <- list(m = m, M_W = M_W)
    of2 <- offs$off2
  } else{
    if(m_memb == T & is.na(sum(M_W)) == F){
      m <- nrow(M_W)
      of2 <- offs$off2
      append <- list(m = m, M_W = M_W)
    }else{
      m <- adj$n ; M_W <- diag(1, adj$n)
      of2 <- offs$off2[1:adj$n]
      append <- list(m = NA, M_W = M_W)
    }
  }
  
  # Covariates
  if(is.null(pars$cov)){
    if(is.null(pars$beta1)) pars$beta2 <- c(1, 1) ; pars$beta1 <- c(1, 1)
    cov <- matrix(0, ncol = length(pars$beta2), nrow = adj$n)
  } else {cov <- as.matrix(pars$cov)}
  
  # Compute relative risks
  rho2 <- M_W%*%(phi_2 + cov%*%pars$beta2 + rep(pars$int[2], adj$n))
  rho1 <- phi_1 + cov%*%pars$beta1 + rep(pars$int[1], adj$n)
  y2 <- rnbinom(m, mu = of2*exp(rho2), size = pars$overd2)
  y1 <- rnbinom(adj$n, mu = offs$off1*exp(rho1), size = pars$overd1)
  
  data_sim <- list(
    # Data
    y1 = y1, y2 = y2, 
    log_offset1 = log(offs$off1), log_offset2 = log(of2),
    X = as.data.frame(pars$cov),
    # Dimensions and adjacency
    W = adj$W, W_n = sum(adj$W)/2,
    k = ncol(as.data.frame(pars$cov)),
    n = adj$n, # m is already in "append"
    # Parameters
    alpha1 = pars$alpha1, alpha2 = pars$alpha2,
    tau1 = pars$tau1, tau2 = pars$tau2,
    phi_1 = phi_1, phi_2 = phi_2,
    mat_phi1 = mat_phi1, mat_phi2 = mat_phi2
  )
  
  return(c(data_sim, append))
}

# Double MM GMCAR-MM
sim_data_Gmcar_double_mm <- function(
  pars, adj, offs, m_memb = F, mcar = F, M_W = NA
){
  library(MASS) # Simulate multivariate normal
  
  # browser()
  # Parameter names
  ## alpha1 - alpha2 - eta0 - eta1 - tau1 - tau2 - 
  ## cov - beta1[2] - beta2[2] - int[2] - overd1 - 
  ## overd 2
  
  # Revert to MCAR from GMCAR
  if(mcar == T){pars$eta1 <- 0 ; pars$alpha2 <- pars$alpha1}
  
  # Covariance matrix (see Jin et. al 2005 paper)
  mat_reg <- pars$eta0*diag(1, adj$n) + pars$eta1*adj$W
  mat_phi1 <- solve(pars$tau1*(adj$D-pars$alpha1*adj$W))
  mat_phi2 <- solve(pars$tau2*(adj$D-pars$alpha2*adj$W))
  sigma_11 <- mat_phi1 + mat_reg%*%mat_phi2%*%mat_reg
  sigma_12 <- mat_reg%*%mat_phi2
  sigma_21 <- mat_phi2%*%mat_reg
  sigma <- rbind(cbind(sigma_11, sigma_12),
                 cbind(sigma_21, mat_phi2))
  
  # Spatial effects
  phi_tot <- mvrnorm(n = 1,
                     mu = c(rep(0, adj$n),
                            rep(0, adj$n)), 
                     Sigma = sigma)
  phi_1 <- phi_tot[1:adj$n]
  phi_2 <- phi_tot[(adj$n + 1):(2*adj$n)]
  
  # Multiple membership
  if(m_memb == T & is.na(sum(M_W)) == T){
    m <- pars$m_memb # number of misaligned areas
    M_W1 <- matrix(runif(adj$n*m), ncol = adj$n, nrow = m) # weights
    M_W1 <- M_W1*matrix(sample(0:1, adj$n*m, replace = T),
                        ncol = adj$n)
    M_W <- M_W1/rowSums(M_W1) # normalized weights
    append <- list(m = m, M_W = M_W)
    of2 <- offs$off2
  } else{
    if(m_memb == T & is.na(sum(M_W)) == F){
      m <- nrow(M_W)
      of2 <- offs$off2
      append <- list(m = m, M_W = M_W)
    }else{
      m <- adj$n ; M_W <- diag(1, adj$n)
      of2 <- offs$off2[1:adj$n]
      append <- list(m = NA, M_W = M_W)
    }
  }
  
  # Covariates
  if(is.null(pars$cov)){
    if(is.null(pars$beta1)) pars$beta2 <- c(1, 1) ; pars$beta1 <- c(1, 1)
    cov <- matrix(0, ncol = length(pars$beta2), nrow = adj$n)
  } else {cov <- as.matrix(pars$cov)}
  
  # Compute relative risks
  rho2 <- M_W%*%(phi_2 + cov%*%pars$beta2) + rep(pars$int[2], m)
  rho1 <- M_W%*%(phi_1 + cov%*%pars$beta1) + rep(pars$int[1], m)
  y2 <- rnbinom(m, mu = of2*exp(rho2), size = pars$overd2)
  y1 <- rnbinom(m, mu = offs$off1*exp(rho1), size = pars$overd1)
  
  data_sim <- list(
    # Data
    y1 = y1, y2 = y2, 
    log_offset1 = log(offs$off1), log_offset2 = log(of2),
    X = as.data.frame(pars$cov),
    # Dimensions and adjacency
    W = adj$W, W_n = sum(adj$W)/2,
    k = ncol(as.data.frame(pars$cov)),
    n = adj$n, # m is already in "append"
    # Parameters
    alpha1 = pars$alpha1, alpha2 = pars$alpha2,
    tau1 = pars$tau1, tau2 = pars$tau2,
    phi_1 = phi_1, phi_2 = phi_2,
    mat_phi1 = mat_phi1, mat_phi2 = mat_phi2
  )
  
  return(c(data_sim, append))
}

# Phi GMCAR
f_phi_gmcar <- function(
  alpha, tau, eta0, eta1, W, mcar = F
){
  library(MASS) # Simulate multivariate normal
  
  # browser()
  # Parameter names
  ## alpha1 - alpha2 - eta0 - eta1 - tau1 - tau2
  
  n <- nrow(W)
  D <- diag(rowSums(W))
  
  # Revert to MCAR from GMCAR
  if(mcar == T){eta1 <- 0 ; alpha[2] <- alpha[1]}
  
  # Covariance matrix (see Jin et. al 2005 paper)
  mat_reg <- eta0*diag(1, n) + eta1*W
  mat_phi1 <- solve(tau[1]*(D-alpha[1]*W))
  mat_phi2 <- solve(tau[2]*(D-alpha[2]*W))
  sigma_11 <- mat_phi1 + mat_reg%*%mat_phi2%*%mat_reg
  sigma_12 <- mat_reg%*%mat_phi2
  sigma_21 <- mat_phi2%*%mat_reg
  sigma <- rbind(cbind(sigma_11, sigma_12),
                 cbind(sigma_21, mat_phi2))
  
  # Spatial effects
  phi_tot <- mvrnorm(n = 1,
                     mu = c(rep(0, n),
                            rep(0, n)), 
                     Sigma = sigma)
  phi_1 <- phi_tot[1:n]
  phi_2 <- phi_tot[(n + 1):(2*n)]
  
  return(list(
    phi1 = phi_1,
    phi2 = phi_2
  ))
}
#

#### FUNCTION TO RETRIEVE RESULTS ####

# Summary for stan objects
f_summary <- function(x){
  c(mean(x), mcse_sd(x), sd(x), quantile(x, c(.025, .1, .5, .9, .975)), 
    Rhat(x), ess_bulk(x))
}

# Check results overall
f_check_spat <- function(fit, pars, real, W, N = 1E4,
                         y1_real, l_RR1_names, l_RR1_real, log_offsets1,
                         y2_real, l_RR2_names, l_RR2_real, log_offsets2,
                         H, biv = F, mm_out = 0
){
  # browser()
  
  # Compute bind and compute estimates - PARAMETERS
  estimates_pars <- data.frame(
    pars =  pars, 
    real = real, 
    mean = sapply(pars, function(x) mean(as.matrix(fit)[, x])),
    median = sapply(pars, function(x) median(as.matrix(fit)[, x]))
  )
  ## Compute Biases - PARAMETERS
  estimates_pars$bias_mean <- (estimates_pars$real - 
                                 estimates_pars$mean)/estimates_pars$real
  estimates_pars$bias_median <- (estimates_pars$real - 
                                   estimates_pars$median)/estimates_pars$real
  
  # Create output list
  output <- list(estimates_pars = estimates_pars)
  
  # Compute bind and compute estimates - RR
  estimates_RR1 <- data.frame(
    RR_real = exp(l_RR1_real),
    RR_est_mean = sapply(l_RR1_names, 
                         function(x) mean(exp(as.matrix(fit)[, x]))),
    RR_est_median = sapply(l_RR1_names, 
                           function(x) median(exp(as.matrix(fit)[, x])))
  )
  
  # Multiple membership outcome or Areal
  if(1 %in% mm_out){
    # Compute bind and compute estimates - Counts
    estimates_counts1 <- data.frame(
      y_obs = y1_real,
      y_est_mean = (H %*% estimates_RR1$RR_est_mean)*exp(log_offsets1),
      y_est_median = (H %*% estimates_RR1$RR_est_median)*exp(log_offsets1)
    )
    
    # Empty list for autocorrelation
    l_auto1 <- list(res1_mean_moran = NULL, res1_median_moran = NULL)
  } else{
    # Compute bind and compute estimates - Counts
    estimates_counts1 <- data.frame(
      y_est_mean = estimates_counts1$RR_est_mean*exp(log_offsets1),
      y_est_median = estimates_counts1$RR_est_median*exp(log_offsets1)
    )
    
    ## Spatial autocorrelation
    res1_mean_moran <- ppt_moran(z = estimates_counts1$y_res_mean, W = W, N = N)
    res1_median_moran <- ppt_moran(z = estimates_counts1$y_res_median, W = W, N = N)
    
    # Save autocorrelation results
    l_auto1 <- list(
      res1_mean_moran = res1_mean_moran, res1_median_moran = res1_median_moran
    )
  }
  
  ## Residuals
  estimates_counts1$y_res_mean <- estimates_counts1$y_obs -
    estimates_counts1$y_est_mean
  estimates_counts1$y_res_median <- estimates_counts1$y_obs -
    estimates_counts1$y_est_median
  ## Bias
  estimates_RR1$bias_RR_mean <- (estimates_RR1$RR_real -
                                   estimates_RR1$RR_est_mean)/estimates_RR1$RR_real
  estimates_RR1$bias_RR_median <- (estimates_RR1$RR_real -
                                     estimates_RR1$RR_est_median)/estimates_RR1$RR_real
  
  if(biv == T){
    # Compute bind and compute estimates - RR
    estimates_RR2 <- data.frame(
      RR_real = exp(l_RR2_real),
      RR_est_mean = sapply(l_RR2_names, 
                           function(x) mean(exp(as.matrix(fit)[, x]))),
      RR_est_median = sapply(l_RR2_names, 
                             function(x) median(exp(as.matrix(fit)[, x])))
    )
    
    # Multiple membership outcome or Areal
    if(2 %in% mm_out){
      # Compute bind and compute estimates - Counts
      estimates_counts2 <- data.frame(
        y_obs = y2_real,
        y_est_mean = (H %*% estimates_RR2$RR_est_mean)*exp(log_offsets2),
        y_est_median = (H %*% estimates_RR2$RR_est_median)*exp(log_offsets2)
      )
      
      # Empty list for autocorrelation
      l_auto2 <- list(res2_mean_moran = NULL, res2_median_moran = NULL)
    } else{
      # Compute bind and compute estimates - Counts
      estimates_counts2 <- data.frame(
        y_est_mean = estimates_counts2$RR_est_mean*exp(log_offsets2),
        y_est_median = estimates_counts2$RR_est_median*exp(log_offsets2)
      )
      ## Spatial autocorrelation
      res2_mean_moran <- ppt_moran(z = estimates_counts2$y_res_mean, W = W, N = N)
      res2_median_moran <- ppt_moran(z = estimates_counts2$y_res_median, W = W, N = N)
      
      # Save autocorrelation results
      l_auto2 <- list(
        res2_mean_moran = res2_mean_moran, res2_median_moran = res2_median_moran
      )
    }
    
    # Residuals
    estimates_counts2$y_res_mean <- estimates_counts2$y_obs -
      estimates_counts2$y_est_mean
    estimates_counts2$y_res_median <- estimates_counts2$y_obs -
      estimates_counts2$y_est_median
    ## Bias
    estimates_RR2$bias_RR_mean <- (estimates_RR2$RR_real -
                                     estimates_RR2$RR_est_mean)/estimates_RR2$RR_real
    estimates_RR2$bias_RR_median <- (estimates_RR2$RR_real -
                                       estimates_RR2$RR_est_median)/estimates_RR2$RR_real
    
    # Save lists
    out_2 <- list(
      estimates_counts2 = estimates_counts2,
      estimates_RR2 = estimates_RR2
    )
    output <- append(output, l_auto2)
    output <- append(output, out_2)
  }
  
  # Append results
  out_1 <- list(
    estimates_counts1 = estimates_counts1,
    estimates_RR1 = estimates_RR1
  )
  output <- append(output, l_auto1)
  output <- append(output, out_1)
  
  return(output)
}
#

#### VARIOUS MATRICES FUNCTIONS ####

# Square root of matrix
f_mat_sqrt <- function(M){
  M_e <- eigen(M)
  if(isSymmetric(M) == T){
    M_e$vectors %*% diag(sqrt(M_e$values)) %*% t(M_e$vectors)
  } else{
    M_e$vectors %*% diag(sqrt(M_e$values)) %*% solve(M_e$vectors)
  }
}

# Find Congruent H
f_cong_H <- function(C_car, C_mm){
  C_sq <- f_mat_sqrt(C_mm %*% C_car)
  C_sq %*% solve(C_car)
}

# Compare matrices with Frobenius norm
f_comp_fr <- function(A, C){
  norm(A - C, type = "F")/norm(C, type = "F")
}
#

#### AUTOCORRELATION TESTS ####

# Geary
f_geary_c <- function(phi, W){
  c_v <- matrix(0, ncol = ncol(W), nrow = nrow(W))
  for(i in 2:nrow(W)){
    for(j in 1:(i-1)){
      c_v[i, j] <- W[i,j]*(phi[i] - phi[j])^2
    }
  }
  (nrow(W) - 1)/(sum(W))*sum(c_v)/sum((phi - mean(phi))^2)
}

# Moran I
moran_test <- function(x, W) {
  n <- length(x)
  z <- as.vector((x - mean(x)))
  (n/sum(W))*(t(z) %*% W %*% z)/sum(z^2)
}
## Permutation
moran_t_perm <- function(x, W, nsim = 1E3) {
  # Size of permutation
  size_perm <- (nrow(W)^2-nrow(W))/2
  W_loop <- W
  
  stat_perm <- numeric()
  
  for(i in 1:nsim){
    W_loop[upper.tri(W_loop)] <- sample(W[upper.tri(W)], size = size_perm)
    W_loop[lower.tri(W_loop)] <- 0
    W_loop <- W_loop + t(W_loop)
    stat_perm[i] <- moran_test(x, W = W_loop)
  }
  
  statistic <- moran_test(x, W)
  
  list(
    statistic = statistic,
    pval = sum(stat_perm > rep(statistic, nsim))/nsim
  )
}
## Permutation test (including a plot).
# ppt_moran <- function(z, W, N=1e4, ...) {
#   stat <- moran(z, W)
#   sim <- replicate(N, moran(sample(z, length(z)), W))
#   p.value <- mean((all <- c(stat, sim)) >= stat)
#   # hist(sim, sub=paste("p =", round(p.value, 4)), xlim=range(all), ...)
#   # abline(v = stat, col="#903030", lty=3, lwd=2)
#   return(p.value)
# }

## Comparsion of the two tests moran (from spdep) and the other
# t_nb <- nb2listw(poly2nb(out_ldn_shp), style = "W")
# yy <- rpois(length(t_nb$weights), lambda = 10)
# 
# col.W <- nb2listw(COL.nb, style="W")
# moran(crime, col.W, length(COL.nb), Szero(col.W))
# moran(yy, listw = t_nb, n = length(t_nb$weights), 
#       S0 =  Szero(t_nb))
# 
# WW <- matrix(0, nrow = length(t_nb$neighbours), 
#              ncol = length(t_nb$neighbours))
# for(i in 1:length(t_nb$neighbours)){
#   pos <- t_nb$neighbours[[i]]
#   wei <- t_nb$weights[[i]]
#   WW[i, pos] <- wei
# }
# 
# moran_test(x = yy, weights = WW)
# 
# ppt_moran(z = yy, W = t_nb)

f_plot_moran <- function(W, y){
  plot(as.numeric(W%*%(y - mean(y))), y)
  abline(lm(y ~ as.numeric(W%*%(y - mean(y)))))
}
#
