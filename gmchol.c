#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void gmchol(double **A, int n, double **R) {
    double **E = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        E[i] = (double *)calloc(n, sizeof(double));
    }
    
    double norm_A = 0.0;
    for (int j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < n; i++) {
            col_sum += fabs(A[i][j]);
        }
        if (col_sum > norm_A) {
            norm_A = col_sum;
        }
    }

    double gamm = 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(A[i][i]) > gamm) {
            gamm = fabs(A[i][i]);
        }
    }

    double delta = fmax(DBL_EPSILON * norm_A, DBL_EPSILON);
    
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
                if (fabs(A[i][j]) > xi_j) {
                    xi_j = fabs(A[i][j]);
                }
            }
        } else {
            xi_j = fabs(A[n - 1][j]);
        }

        double beta_j = sqrt(fmax(gamm, fmax(xi_j / n, DBL_EPSILON)));

        if (delta >= fmax(fabs(phi_j), (theta_j * theta_j) / (beta_j * beta_j))) {
            E[j][j] = delta - phi_j;
        } else if (fabs(phi_j) >= fmax((delta * delta) / (beta_j * beta_j), delta)) {
            E[j][j] = fabs(phi_j) - phi_j;
        } else if ((theta_j * theta_j) / (beta_j * beta_j) >= fmax(delta, fabs(phi_j))) {
            E[j][j] = ((theta_j * theta_j) / (beta_j * beta_j)) - phi_j;
        }

        R[j][j] = sqrt(A[j][j] - sum_val + E[j][j]);
    }

    // Free allocated memory for E
    for (int i = 0; i < n; i++) {
        free(E[i]);
    }
    free(E);
}

// Example
int main() {
    int n = 3; // Matrix size
    double **A = (double **)malloc(n * sizeof(double *));
    double **R = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
        R[i] = (double *)calloc(n, sizeof(double)); // Init R to zero
    }

    // Symmetric positive definite matrix
    A[0][0] = 4; A[0][1] = 12; A[0][2] = -16;
    A[1][0] = 12; A[1][1] = 37; A[1][2] = -43;
    A[2][0] = -16; A[2][1] = -43; A[2][2] = 98;

    gmchol(A, n, R);

    printf("R matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", R[i][j]);
        }
        printf("\n");
    }

    // Free memory allocated for A and R
    for (int i = 0; i < n; i++) {
        free(A[i]);
        free(R[i]);
    }
    free(A);
    free(R);

    return 0;
}

