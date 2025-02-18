# R port of the Gill/Murray cholesky Decomposition
# Input matrix A must be non-singular and symmetric
# Author: Wyatt Gill. Credits to Jeff Gill and Gary King

gmchol <- function(A) {

  n <- nrow(A)
  R <- diag(n)
  E <- matrix(0, n, n)
  norm_A <- max(colSums(abs(A)))
  gamm <- max(abs(diag(A)))
  delta <- max(.Machine$double.eps * norm_A, .Machine$double.eps)
  
  for (j in 1:n) {
    theta_j <- 0
    
    for (i in 1:n) {
      sum_val <- sum(R[1:(i-1), i] * R[1:(i-1), j], na.rm = TRUE)
      R[i, j] <- (A[i, j] - sum_val) / R[i, i]
      
      if ((A[i, j] - sum_val) > theta_j) {
        theta_j <- A[i, j] - sum_val
      }
      
      if (i > j) {
        R[i, j] <- 0
      }
    }
    
    sum_val <- sum(R[1:(j-1), j]^2, na.rm = TRUE)
    phi_j <- A[j, j] - sum_val
    
    if ((j + 1) <= n) {
      xi_j <- max(abs(A[(j+1):n, j]), na.rm = TRUE)
    } else {
      xi_j <- abs(A[n, j])
    }
    
    beta_j <- sqrt(max(gamm, xi_j / n, .Machine$double.eps))
    
    if (delta >= max(abs(phi_j), (theta_j^2) / (beta_j^2))) {
      E[j, j] <- delta - phi_j
    } else if (abs(phi_j) >= max((delta^2) / (beta_j^2), delta)) {
      E[j, j] <- abs(phi_j) - phi_j
    } else if ((theta_j^2) / (beta_j^2) >= max(delta, abs(phi_j))) {
      E[j, j] <- ((theta_j^2) / (beta_j^2)) - phi_j
    }
    
    R[j, j] <- sqrt(A[j, j] - sum_val + E[j, j])
  }
  
  return(t(R) %*% R)
}
