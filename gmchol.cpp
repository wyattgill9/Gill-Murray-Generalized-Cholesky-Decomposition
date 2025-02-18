#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

void gmchol(const std::vector<std::vector<double>>& A, int n, std::vector<std::vector<double>>& R) {
    std::vector<std::vector<double>> E(n, std::vector<double>(n, 0.0));
    
    double norm_A = 0.0;
    for (int j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < n; i++) {
            col_sum += std::fabs(A[i][j]);
        }
        if (col_sum > norm_A) {
            norm_A = col_sum;
        }
    }

    double gamm = 0.0;
    for (int i = 0; i < n; i++) {
        if (std::fabs(A[i][i]) > gamm) {
            gamm = std::fabs(A[i][i]);
        }
    }

    double delta = std::fmax(std::numeric_limits<double>::epsilon() * norm_A, std::numeric_limits<double>::epsilon());
    
    for (int j = 0; j < n; j++) {
        double theta_j = 0.0;

        for (int i = 0; i < n; i++) {
            double sum_val = 0.0;
            for (int k = 0; k < i; k++) {
                sum_val += R[k][i] * R[k][j];
            }
            R[i][j] = (A[i][j] - sum_val) / R[i][i];

            if ((A[i][j] - sum_val) > theta_j) {
                theta_j = A[i][j] - sum_val;
            }

            if (i > j) {
                R[i][j] = 0.0;
            }
        }

        double sum_val = 0.0;
        for (int k = 0; k < j; k++) {
            sum_val += R[k][j] * R[k][j];
        }
        double phi_j = A[j][j] - sum_val;

        double xi_j;
        if ((j + 1) < n) {
            xi_j = 0.0;
            for (int i = j + 1; i < n; i++) {
                if (std::fabs(A[i][j]) > xi_j) {
                    xi_j = std::fabs(A[i][j]);
                }
            }
        } else {
            xi_j = std::fabs(A[n - 1][j]);
        }

        double beta_j = std::sqrt(std::fmax(gamm, std::fmax(xi_j / n, std::numeric_limits<double>::epsilon())));

        if (delta >= std::fmax(std::fabs(phi_j), (theta_j * theta_j) / (beta_j * beta_j))) {
            E[j][j] = delta - phi_j;
        } else if (std::fabs(phi_j) >= std::fmax((delta * delta) / (beta_j * beta_j), delta)) {
            E[j][j] = std::fabs(phi_j) - phi_j;
        } else if ((theta_j * theta_j) / (beta_j * beta_j) >= std::fmax(delta, std::fabs(phi_j))) {
            E[j][j] = ((theta_j * theta_j) / (beta_j * beta_j)) - phi_j;
        }

        R[j][j] = std::sqrt(A[j][j] - sum_val + E[j][j]);
    }
}

// Example
int main() {
    int n = 3; // Matrix size
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0)); // Init R to zero

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

