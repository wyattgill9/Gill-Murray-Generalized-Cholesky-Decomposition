#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void gmchol(double **A, int n, double **R) {
    double norm_A = 0.0, gamm = 0.0;

    for (int j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int i = 0; i < n; i++) {
            double val = fabs(A[i][j]);
            col_sum += val;
            if (i == j && val > gamm) gamm = val;
        }
        if (col_sum > norm_A) norm_A = col_sum;
    }

    double delta = fmax(DBL_EPSILON * norm_A, DBL_EPSILON);

    for (int j = 0; j < n; j++) {
        double theta_j = 0.0;

        for (int i = 0; i <= j; i++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += R[k][i] * R[k][j];
            }

            if (i < j) {
                R[i][j] = (A[i][j] - sum) / R[i][i];
                double residual = A[i][j] - sum;
                if (fabs(residual) > theta_j) theta_j = fabs(residual);
            } else {
                double phi_j = A[j][j];
                for (int k = 0; k < j; k++) {
                    phi_j -= R[k][j] * R[k][j];
                }

                double xi_j = 0.0;
                for (int t = j + 1; t < n; t++) {
                    double val = fabs(A[t][j]);
                    if (val > xi_j) xi_j = val;
                }

                double beta_j = sqrt(fmax(gamm, fmax(xi_j / n, DBL_EPSILON)));
                double threshold = fmax(fabs(phi_j), (theta_j * theta_j) / (beta_j * beta_j));

                double correction = 0.0;
                if (delta >= threshold) {
                    correction = delta - phi_j;
                } else if (fabs(phi_j) >= fmax((delta * delta) / (beta_j * beta_j), delta)) {
                    correction = fabs(phi_j) - phi_j;
                } else {
                    correction = ((theta_j * theta_j) / (beta_j * beta_j)) - phi_j;
                }

                R[j][j] = sqrt(phi_j + correction);
            }
        }

        for (int i = j + 1; i < n; i++) {
            R[i][j] = 0.0;
        }
    }
}

int main() {
    int n = 3;
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

