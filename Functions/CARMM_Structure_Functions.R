#### CONDITIONING FROM COVARIANCE FUNCTION ####

# Ver Hoef
## Recover structure - matrix R
calcMatrixR <- function(matrixSigma, matrixH){
  
  matrixSigmaTilde <- (matrixH %*% matrixSigma) %*% t(matrixH)
  matrixQTilde <- solve(matrixSigmaTilde)
  matrixQTilde <- ginv(matrixSigmaTilde)
  matrixS <- diag(diag(matrixQTilde))
  matrixR <- matrixS-matrixQTilde
  
}
## Recover structure - precision matrix
calcMatrixQtilde <- function(matrixSigma, matrixH){
  library(MASS)
  
  matrixSigmaTilde <- (matrixH %*% matrixSigma) %*% t(matrixH)
  matrixQTilde <- ginv(matrixSigmaTilde)
  
  return(matrixQTilde)
}

# Gaussian conditionals (USES MOORE-PENROSE INVERSE)
f_cond_pair <- function(C_Q){
  library(MASS)
  
  # browser()
  
  # Base lists
  ## Sub-Covariance (Sigma_XX)
  S_XX <- list()
  ## Cross sub-covariance S_XY
  S_XY <- list()
  ## Conditioning covariance (Sigma_YY)
  S_YY <- list()
  ## Number of variables
  n <- nrow(C_Q)
  
  # OFF DIAGONAL - CONDITIONING QUANTITIES
  for(i in 2:n){
    # Base lists for loop
    s_t_xx <- list()
    s_t_xy <- list()
    s_t_yy <- list()
    for(j in 1:(i-1)){
      # Conditioning quantities
      s_t_xx[[j]] <- C_Q[c(i, j), c(i, j)]
      s_t_xy[[j]] <- C_Q[c(i, j), -c(i, j)]
      s_t_yy[[j]] <- C_Q[-c(i, j), -c(i, j)]
    }
    # Save lists
    S_XX[[i]] <- s_t_xx
    S_XY[[i]] <- s_t_xy
    S_YY[[i]] <- s_t_yy
  }
  
  # CONDITIONAL COVARINACE AND VARIANCE
  ## Covariance
  C_XY <- matrix(0, ncol = n, nrow = n)
  ## Variance
  V_XY <- matrix(0, ncol = n, nrow = n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      
      s_t_cxy <- S_XX[[i]][[j]] -
        S_XY[[i]][[j]]%*%ginv(S_YY[[i]][[j]])%*%t(S_XY[[i]][[j]])
      
      # Results
      ## Covariance
      C_XY[i, j] <- s_t_cxy[1, 2]
      ## Variance
      V_XY[i, j] <- s_t_cxy[1, 1]
      V_XY[j, i] <- s_t_cxy[2, 2]
    }
  }
  
  # FINAL CORRELATION
  R_XY <- matrix(0, ncol = n, nrow = n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      if(V_XY[i, j] <= 0 | V_XY[j, i] <= 0){
        R_XY[i, j] <- 0
      } else{
        R_XY[i, j] <- C_XY[i, j]/sqrt(V_XY[i, j]*V_XY[j, i])
      }
    }
  }
  
  # Impose symmetry
  ## COVARIANCE
  C_XY <- C_XY + t(C_XY)
  ## VARIANCE
  # V_XY <- V_XY + t(V_XY)
  ## CORRELATION
  R_XY <- R_XY + t(R_XY)
  
  # Diagonal - CORRELATION
  diag(R_XY) <- rep(1, n)
  
  output <- list(
    C_XY = C_XY, R_XY = R_XY, V_XY = V_XY
  )
  return(output)
}

# Dummify matrix 
f_mat_dummy <- function(M, d){
  M <- round(M, d)
  M[which(M != 0)] <- 1
  return(M)
}

# FUll membership CARMM specification
f_carmm_sp_prec <- function(W, H, alpha, tau){
  K <- solve(H)
  D <- diag(rowSums(W))
  
  D_t <- t(K) %*% D %*% K
  W_t <- t(K) %*% W %*% K
  
  tau*(D_t - alpha*W_t)
}
#

#### CORRELATION FUNCTIONS ####

# Fisher z-transform
## page 149 PenPC: A Two-Step Approach to Estimate the Skeletons of
## High-Dimensional Directed Acyclic Graphs
z_fisher <- function(rho) .5*log((1+rho)/(1-rho))

# Partial correlation - Individual
f_parCor <- function(x, y, Z){
  # Control necessary for single neighbour area
  if(length(Z) != 0){
    r_x <- residuals(lm(x ~ Z))
    r_y <- residuals(lm(y ~ Z))
    
    # Return correlation
    cor(r_x, r_y)
  } else{
    cor(x, y)
  }
}

# Partial correlation - Aggregate
car_parCor <- function(phi, W, complete = F){
  # Number of areas
  K <- nrow(W)
  
  # Number of samples
  n <- nrow(phi)
  
  # Correlation matrix
  phi_cor <- matrix(NA, ncol = K, nrow = K)
  
  # Test statistics
  phi_test <- matrix(NA, ncol = K, nrow = K)
  
  # Conditional correlation
  if(complete == T){
    # Squared normalising factor for test (global since all areas considered)
    t_factor <- n - (K-2) - 3
    
    # Base matrices
    phi_cor <- matrix(0, nrow = K, ncol = K)
    phi_test <- phi_cor
    
    for(i in 2:K){
      for(j in 1:(i-1)){
        
        # Partial correlation
        phi_cor[i, j] <- f_parCor(phi[,i], phi[,j], phi[, -c(i,j)])
        
        # Test statistic
        phi_test[i, j] <- abs(z_fisher(phi_cor[i, j]))*sqrt(t_factor)
      }
    }
    
    # Symmetric matrices
    ## Cor
    phi_cor <- phi_cor + t(phi_cor)
    diag(phi_cor) <- rep(1, K)
    ## Test
    phi_test <- phi_test + t(phi_test)
    
    
  } else{
    for(i in 1:K){
      # Squared normalising factor for test
      t_factor <- n - length(phi[i,which(W[i,] == 1)]) - 3
      
      # Non-neighbours
      phi_cor[i, -c(which(W[i,] == 1), i)] <-
        apply(phi[, -c(which(W[i,] == 1), i)], 2, 
              function(x) f_parCor(phi[,i], x, phi[,which(W[i,] == 1)])
        )
      ## Test statistic
      phi_test[i, -c(which(W[i,] == 1), i)] <-
        abs(z_fisher(phi_cor[i, -c(which(W[i,] == 1), i)]))*sqrt(t_factor)
      
      # Neighbours
      for(j in which(W[i,] == 1)){
        phi_cor[i, j] <-
          f_parCor(phi[,i], phi[, j], phi[, -c(j, which(W[i,] == 0))])
        phi_test[i, j] <- abs(z_fisher(phi_cor[i, j]))*sqrt(t_factor - 1)
      }
      
      # Diagonal
      phi_cor[i, i] <- 1
      phi_test[i, i] <- 0
    }
  }
  
  output <- list(phi_cor = phi_cor, phi_test = phi_test)
  return(output)
}

# Estimated adj matrix with partial correlations
est_carStr <- function(phi_corMat, threshold, W){
  rawLook <- phi_corMat
  rawLook[which(abs(rawLook) < threshold)] <- 0
  diag(rawLook) <- 0
  rawLook[which(abs(rawLook) >= threshold)] <- 1
  
  output <- list(mat_estStr = rawLook, n_err = sum(W - rawLook))
  return(output)
}

# LASSO partial correlation
f_lasso_str <- function(W, phi, check = F){
  library(glmnet)
  
  # Base vectors
  n_missed <- numeric()
  n_extra <- numeric()
  non_z_l <- list()
  
  # Change names
  colnames(phi) <- 1:ncol(phi)
  
  if(check == F){
    for(i in 1:ncol(phi)){
      # browser()
      cv_phi <- cv.glmnet(x = phi[, -i], y = phi[, i], alpha = 1)
      
      lasso_coef <- predict(cv_phi, type = "coefficients")
      
      non_z_l[[i]] <- as.numeric(names(which(lasso_coef[,1] != 0))[-1])
    }
    
    # Return
    output <- list(lasso_str = non_z_l)
    
  } else{
    for(i in 1:ncol(phi)){
      cv_phi <- cv.glmnet(x = phi[, -i], y = phi[, i], alpha = 1)
      
      lasso_coef <- predict(cv_phi, type = "coefficients")
      
      non_z_l[[i]] <- as.numeric(names(which(lasso_coef[,1] != 0))[-1])
      
      # Errors
      n_missed[i] <- sum(!which(W[i, ] == 1)%in%non_z_l[[i]])
      n_extra[i] <- sum(!non_z_l[[i]]%in%which(W[i, ] == 1))
    }
    
    # Return
    output <- list(
      lasso_str = non_z_l,
      n_missed = n_missed,
      n_extra = n_extra
    )
  }
  
  return(output)
}

# Compare two adjacency matrices
f_comp_mat <- function(W1, W2){
  # Comparison matrix
  W <- W1*2 - W2
  
  ggplot(data.frame(x = rep(1:nrow(W), each = nrow(W)),
                    y = rep(1:nrow(W), nrow(W)),
                    z = as.numeric(W)),
         aes(x, y, fill = as.character(z))) +
    geom_tile() + scale_y_reverse()
}

# Correlation matrix from known covariance
f_cor_mat <- function(S){
  D_s_i <- diag(1/sqrt(diag(S)))
  D_s_i %*% S %*% D_s_i
}
#

#### PARTITION MATRIX ####

# Generate partitioned matrix
f_part_mat_2 <- function(A, r = 12){
  library(Matrix)
  
  # Checks
  ## Dimension control
  if(nrow(A) != ncol(A)){
    stop("Rectangular matrix")
  }
  ## Full rank control
  if(rankMatrix(A)[1] == nrow(A)){
    stop("Full rank matrix")
  }
  
  # Find linear dependent columns
  ## QR - R part
  R_A <- round(qr.R(qr(A, tol = 1E-19)), r)
  ## Linear dependent columns
  ld_cols <- which(diag(R_A) == 0)
  
  # Full rank partition
  A_11 <- A[-ld_cols, -ld_cols]
  
  # Compute B
  A_1 <- A[, -ld_cols]
  A_2 <- A[, ld_cols]
  Q_1 <- qr.Q(qr(A_1, tol = 1E-19), complete = T)[,-ld_cols]
  R_1 <- qr.R(qr(A_1, tol = 1E-19), complete = T)[-ld_cols,]
  B <- solve(R_1) %*% t(Q_1) %*% A_2
  
  # Remaining partition
  A_off_diag <- A_11 %*% B
  A_diag <- t(B) %*% A_off_diag
  
  # Final matrix
  final_mat <- rbind(cbind(A_11, A_off_diag), cbind(t(A_off_diag), A_diag))
  
  # Rank
  r_final <- rankMatrix(final_mat)[1]
  
  # General control
  if(rankMatrix(A) != r_final){
    stop("Initial and final matrices ranks differ")
  }
  
  # Output
  output <- list(
    full = final_mat,
    rank = r_final,
    non_sing = final_mat[1:r_final, 1:r_final],
    B_mat = B
  )
}

# Second aspect
f_part_mat <- function(A, r = 12){
  library(Matrix)
  
  browser()
  
  # Checks
  ## Dimension control
  if(nrow(A) != ncol(A)){
    stop("Rectangular matrix")
  }
  ## Full rank control
  rank_A <- rankMatrix(A)[1]
  if(rank_A == nrow(A)){
    stop("Full rank matrix")
  }
  
  # Find linear dependent columns
  ld_cols <- numeric()
  j <- 1
  for(i in 1:nrow(A)){
    if(rankMatrix(A[, -i])[1] != rank_A){
      ld_cols[j] <- i
      j <- j + 1
    }
  }
  
  # Full rank partition
  A_11 <- A[-ld_cols, -ld_cols]
  
  # Compute B
  A_1 <- A[, -ld_cols]
  A_2 <- A[, ld_cols]
  Q_1 <- qr.Q(qr(A_1, tol = 1E-19), complete = T)[,-ld_cols]
  R_1 <- qr.R(qr(A_1, tol = 1E-19), complete = T)[-ld_cols,]
  B <- solve(R_1) %*% t(Q_1) %*% A_2
  
  # Remaining partition
  A_off_diag <- A_11 %*% B
  A_diag <- t(B) %*% A_off_diag
  
  # Final matrix
  final_mat <- rbind(cbind(A_11, A_off_diag), cbind(t(A_off_diag), A_diag))
  
  # Rank
  r_final <- rankMatrix(final_mat)[1]
  
  # General control
  if(rankMatrix(A) != r_final){
    stop("Initial and final matrices ranks differ")
  }
  
  # Output
  output <- list(
    full = final_mat,
    rank = r_final,
    non_sing = final_mat[1:r_final, 1:r_final],
    B_mat = B
  )
}
#

#### SYMMETRIC SQUARE MATRICES HEATMAP FUNCTION ####

## NOTE ON ORDERING ELEMENTS IN HEATMAPS AND MATRICES
### expand.grid : orders left-right by row
### as.numeric(matrix) : left-right by column
### Positionally the two methods match, the problem is the matrix coords
### arrangement. So, in the plot "row" has to be the the y-axis and "column" 
### the x-axis. In order to match plot and matrix coords we need to use
### as.numeric(t(c_mat)), so that it's printed by row instead of column.
### The 1st argument in expand.grid has to be the rows and the 2nd the columns
## Example
# c_mat <- matrix(1:9, ncol = 3, nrow = 3) ; as.numeric(c_mat)
# c_grid <- expand.grid(row = 1:3, column = 1:3)
# c_grid$vals <- as.numeric(t(c_mat)) # now they match in terms of the mat coords
# ggplot(c_grid, aes(x = column, y = row, fill = as.factor(vals))) +
#   scale_y_reverse() + geom_tile()

# Function for heatmatrix
## ATTENTION: ONLY TESTED W/ SYMMETRIC SQUARE MATRICES!!!
## SHOULD ALSO WORK FOR CAR GENERATED ON A GRID
heat_matrix <- function(M){ 
  library(ggplot2)
  
  # Create grid
  side <- nrow(M)
  base_grid <- expand.grid(row = 1:side, column = 1:side)
  base_grid$vals <- as.numeric(t(M)) # now they match in terms of matr coords
  
  # Plot
  ggplot(base_grid, aes(x = column, y = row, fill = vals)) +
    scale_y_reverse() + geom_tile()
}
#

#### DIRECT PARTIAL CORRELATION FUNCTION ####

# Partial correlation with generalised inverse
## similar to library(corpcor) ; cor2pcor(cov2cor(H %*% S %*% t(H)))
partial_cor_ginv <- function(S, H, zero_threshold = 0){
  # browser()
  S_mm <- H %*% S %*% t(H)
  S_inv_mm <- ginv(S_mm)
  D_p_inv <- diag(1/sqrt(diag(S_inv_mm)))
  R_p <- - D_p_inv %*% S_inv_mm %*% D_p_inv
  diag(R_p) <- rep(1, nrow(R_p))
  
  R_p_inf <- adj_inferred(abs(R_p) > zero_threshold)
  
  return(
    list(R_p = R_p, R_p_inf = R_p_inf)
  )
}

# Direct partial correlation function
# f_pcor <- function(m){
#   # browser()
#   
#   library(MASS)
#   
#   # Compute CORRELATION matrix
#   r <- diag(sqrt(1/diag(m))) %*% m %*% diag(sqrt(1/diag(m)))
# 
#   # Compute PARTIAL CORRELATION matrix
#   r_inv <- ginv(r)
#    diag(sqrt(1/diag(r_inv))) %*% r_inv %*% diag(sqrt(1/diag(r_inv)))
#   
#   # m <- -ginv(m)
#   # diag(m) <- -diag(m)
#   # m
# }
#
