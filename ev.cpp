#include <stdio.h>
#include "functions.h"
#include <math.h>
#include "stdlib.h"


double norm_matrix(int n, double *A) {
    double norm = 0, sum = 0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            sum += fabs(A[j * n + i]);
        }
        if(sum > norm) {
            norm = sum;
        }
        sum = 0;
    }
    return norm;
}


double Nevyazka_1(int n, double *A, double *X) {



    double A_sum = 0, X_sum = 0;
    for(int i = 0; i < n; i++) {
        A_sum += A[i * n + i];
        X_sum += X[i];
    }
    return fabs(A_sum - X_sum);
}


double Nevyazka_2(int n, double *A, double *X) {

    double sum = 0, X_sum = 0;
    for(int i = 0; i < n; i++) {
        X_sum += X[i] * X[i];
        for(int j = 0; j < n; j++) {
            sum += A[i * n + j] * A[i * n + j];
        }
    }
    return fabs(sqrt(sum) - sqrt(X_sum));
}

void pochti_triangle(int n, double *A, double *x) {
    double EPS = 1e-100;
    double tmp1 = 0;
    double tmp2 = 0;
    for(int i = 0; i < n - 2; i++) {
        tmp1 = 0;
        for(int j = i + 2; j < n; j++) {
            tmp1 += A[j * n + i] * A[j * n + i];
        }

        tmp2 = sqrt(A[(i + 1) * n + i] * A[(i + 1) * n + i] + tmp1);

        if(tmp2 < EPS) {
            A[(i + 1) * n + i] = 0;
            A[(i + 2) * n + i] = 0;
            continue;
        }
        if(tmp1 < EPS) {
            A[(i + 2) * n + i] = 0;
            continue;
        }

        x[i + 1] = A[(i + 1) * n + i] - tmp2;
        tmp1 = 1 / sqrt(x[i + 1] * x[i + 1] + tmp1);
        x[i + 1] *= tmp1;

        for(int j = i + 2; j < n; j++) {
            x[j] = A[j * n + i] * tmp1;
        }

        for(int j = i + 1; j < n; j++) {
            tmp1 = 0;
            for(int k = i + 1; k < n; k++) {
                tmp1 += x[k] * A[k * n + j];
            }
            tmp1 *= 2;
            for(int k = i + 1; k < n; k++)
                A[k * n + j] -= tmp1 * x[k];
        }
        for(int j = 0; j < n; j++) {
            tmp1 = 0;
            for(int k = i + 1; k < n; k++) {
                tmp1 += x[k] * A[j * n + k];
            }
            tmp1 *= 2;
            for(int k = i + 1; k < n; k++) {
                A[j * n + k] -= tmp1 * x[k];
            }
        }
        A[(i + 1) * n + i] = tmp2;

    }
}



double dvigai_value(int n, double *A, int k) {
    return A[(k - 1) * n + k - 1] + 0.5 * A[(k - 1) * n + k - 2];
    //return A[(k - 1) * n + k - 1];
}


void dvigai(int n, double *A, int k, double s) {
    for(int i = 0; i < k; i++) {
        A[i * n + i] -= s;
    }
}

void dvigai_back(int n, double *A, int k, double s) {
    for(int i = 0; i < k; i++) {
        A[i * n + i] += s;
    }
}

void do_QR(int n, double *A, int k, double *x1, double *x2) {
    double EPS = 1e-100;
    double tmp1 = 0;
    double tmp2 = 0;
    for(int i = 0; i < k - 1; i++) {
        tmp1 = A[(i + 1) * n + i] * A[(i + 1) * n + i];
        if(tmp1 < EPS) {
            x1[i] = 0;
            x2[i] = 0;
            continue;
        }

        tmp2 = sqrt(A[i * n + i] * A[i * n + i] + tmp1) * ((A[i * n + i] > 0) ? -1 : 1);
        A[i * n + i] -= tmp2;
        tmp1 = sqrt(A[i * n + i] * A[i * n + i] + tmp1);
        x1[i] = A[i * n + i] / tmp1;
        x2[i] = A[(i + 1) * n + i] / tmp1;

        for(int j = i + 1; j < k; j++) {
            tmp1 = x1[i] * A[i * n + j];
            tmp1 += x2[i] * A[(i + 1) * n + j];

            tmp1 *= 2;

            A[i * n + j] -= tmp1 * x1[i];
            A[(i + 1) * n + j] -= tmp1 * x2[i];
        }
        A[i * n + i] = tmp2;
        A[(i + 1) * n + i] = 0;
    }
}


void do_RQ(int n, double *A, int k, double *x1, double *x2) {
    double EPS = 1e-100;
    double tmp = 0;
    for(int i = 0; i < k - 1; i++) {
        if (fabs (x1[i]) < EPS && fabs (x2[i]) < EPS) {
            continue;
        }

        for(int j = 0; j < i + 2; j++) {
            tmp = A[j * n + i] * x1[i];
            tmp += A[j * n + i + 1] * x2[i];
            tmp *= 2;
            A[j * n + i] -= tmp * x1[i];
            A[j * n + i + 1] -= tmp * x2[i];
        }
    }
}

int find_sobsn_znach(int n, double *A, double *X, double eps, double *x1, double *x2) {
    int ev_iter = 0;
    double tmp = 0;
    double s = 0;

    for(int k = n; k > 2; k--) {
        ev_iter = 0;

        while(fabs(A[(k - 1) * n + k - 2]) > eps) {
            s = dvigai_value(n, A, k);
            dvigai(n, A, k, s);
            do_QR(n, A, k, x1, x2);

            do_RQ(n, A, k, x1, x2);
            dvigai_back(n, A, k, s);

            ev_iter++;
            if(ev_iter > 100) {

                return 1;
            }
        }
    }

    if(n > 1) {
        tmp = A[0 * n + 0] + A[1 * n + 1];
        s = A[0 * n + 0] * A[1 * n + 1] - A[0 * n + 1] * A[1 * n + 0];
        s = tmp * tmp - 4.0 * s;
        if (s < 0.0) {

            return 2;
        }

        s = sqrt (s);

        A[0 * n + 0] = 0.5 * (tmp + s);
        A[1 * n + 1] = 0.5 * (tmp - s);
    }

    for(int i = 0; i < n; ++i)
        X[i] = A[i * n + i];


    return 0;
}
