#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

/**
 * Gill-Murray Modified Cholesky Decomposition
 * 
 * Reference: Gill, Jeff and Gary King. "What to do When Your Hessian is Not Invertible:
 * Alternatives to Model Respecification in Nonlinear Estimation," Sociological
 * Methods and Research, Vol. 32, No. 1 (2004): Pp. 54-87.
 * 
 * Input: A - symmetric matrix (n x n)
 * Output: Modified matrix A + E where E is a diagonal correction matrix
 * Returns: R'R where R is upper triangular
 */
std::vector<std::vector<double>> gmchol(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        R[i][i] = 1.0;
    }
    
    std::vector<std::vector<double>> E(n, std::vector<double>(n, 0.0));
    
    double norm_A = 0.0;
    for (int j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < n; i++) {
            col_sum += std::fabs(A[i][j]);
        }
        norm_A = std::max(norm_A, col_sum);
    }
    
    double gamm = 0.0;
    for (int i = 0; i < n; i++) {
        gamm = std::max(gamm, std::fabs(A[i][i]));
    }
    
    double eps = std::numeric_limits<double>::epsilon();
    double delta = std::max(eps * norm_A, eps);
    
    for (int j = 0; j < n; j++) {
        double theta_j = 0.0;
        
        for (int i = 0; i <= j; i++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += R[k][i] * R[k][j];
            }
            
            if (i < j) {
                R[i][j] = (A[i][j] - sum) / R[i][i];
                if ((A[i][j] - sum) > theta_j) {
                    theta_j = A[i][j] - sum;
                }
            }
        }
        
        double sum = 0.0;
        for (int k = 0; k < j; k++) {
            sum += R[k][j] * R[k][j];
        }
        double phi_j = A[j][j] - sum;
        
        double xi_j = 0.0;
        if (j + 1 < n) {
            for (int i = j + 1; i < n; i++) {
                xi_j = std::max(xi_j, std::fabs(A[i][j]));
            }
        } else {
            xi_j = std::fabs(A[n-1][j]);
        }

        double beta_j = std::sqrt(std::max({gamm, xi_j / n, eps}));
        
        double abs_phi_j = std::fabs(phi_j);
        double theta_sq_over_beta_sq = (theta_j * theta_j) / (beta_j * beta_j);
        double delta_sq_over_beta_sq = (delta * delta) / (beta_j * beta_j);
        
        if (delta >= std::max(abs_phi_j, theta_sq_over_beta_sq)) {
            E[j][j] = delta - phi_j;
        } else if (abs_phi_j >= std::max(delta_sq_over_beta_sq, delta)) {
            E[j][j] = abs_phi_j - phi_j;
        } else if (theta_sq_over_beta_sq >= std::max(delta, abs_phi_j)) {
            E[j][j] = theta_sq_over_beta_sq - phi_j;
        }
        
        R[j][j] = std::sqrt(A[j][j] - sum + E[j][j]);
    }
    
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += R[k][i] * R[k][j];
            }
            result[i][j] = sum;
        }
    }
    
    return result;
}

// Example
int main() {
    int n = 3; // Matrix size
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0)); 

    // Symmetric positive definite matrix
    A[0][0] = 4; A[0][1] = 12; A[0][2] = -16;
    A[1][0] = 12; A[1][1] = 37; A[1][2] = -43;
    A[2][0] = -16; A[2][1] = -43; A[2][2] = 98;

    gmchol(A, n, R);

    std::cout << "R matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << R[i][j] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}

