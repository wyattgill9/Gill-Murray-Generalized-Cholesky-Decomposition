# Python port of the Gill/Murray cholesky Decomposition
# Input matrix A must be non-singular and symmetric
# Author: Wyatt Gill. Credits to Jeff Gill and Gary King

import numpy as np

def gmchol(A):
    n = A.shape[0]
    R = np.eye(n)
    E = np.zeros((n, n))
    norm_A = np.max(np.sum(np.abs(A), axis=0))
    gamm = np.max(np.abs(np.diag(A)))
    delta = max(np.finfo(float).eps * norm_A, np.finfo(float).eps)
    
    for j in range(n):
        theta_j = 0
        
        for i in range(n):
            sum_val = np.sum(R[:i, i] * R[:i, j])
            R[i, j] = (A[i, j] - sum_val) / R[i, i]
            
            if (A[i, j] - sum_val) > theta_j:
                theta_j = A[i, j] - sum_val
            
            if i > j:
                R[i, j] = 0
        
        sum_val = np.sum(R[:j, j] ** 2)
        phi_j = A[j, j] - sum_val
        
        if (j + 1) < n:
            xi_j = np.max(np.abs(A[(j + 1):n, j]))
        else:
            xi_j = np.abs(A[n - 1, j])
        
        beta_j = np.sqrt(max(gamm, xi_j / n, np.finfo(float).eps))
        
        if delta >= max(np.abs(phi_j), (theta_j ** 2) / (beta_j ** 2)):
            E[j, j] = delta - phi_j
        elif np.abs(phi_j) >= max((delta ** 2) / (beta_j ** 2), delta):
            E[j, j] = np.abs(phi_j) - phi_j
        elif (theta_j ** 2) / (beta_j ** 2) >= max(delta, np.abs(phi_j)):
            E[j, j] = ((theta_j ** 2) / (beta_j ** 2)) - phi_j
        
        R[j, j] = np.sqrt(A[j, j] - sum_val + E[j, j])
    
    return np.dot(R.T, R)

