#include <stdio.h>
#include <time.h>
#include "stdlib.h"
#include "functions.h"
#include <math.h>
#include <new>
#include <climits>
int main(int argc, char *argv[]) {
    int n = 0, m = 0, s = 0;
    double eps = 0;
    char *name = nullptr;
    double *A = nullptr;
    double r1 = -1, r2 = -1, t2 = 0;
    int fl_1 = 0;

    double *X = nullptr;
    double *x1 = nullptr;
    double *x2 = nullptr;
    FILE *input = nullptr;

    if(!((argc == 5 || argc == 6) && (sscanf(argv[1], "%d", &n) == 1) && (sscanf(argv[2], "%d", &m) == 1) && (sscanf(argv[3], "%lf", &eps) == 1)&& (sscanf(argv[4], "%d", &s) == 1)&& ((s >= 0 && s <= 4)|| (s == 0 && argc == 5)))) {
        printf("Data input error");
        return -1;
    }

    if( n - sqrt(INT_MAX)>= 0)
    {
        printf("n*n > int_max. Cannot allocate!\n");
        return 2;
    }

  A = new double[n * n];
    if(!A) {
        printf("Memory allocation error");
        return -1;
    }
    X = new double[n];
    if(!X) {
	printf("Memory allocation error");
		delete []A;
		return -1;
    }
    x1 = new double[n];
    if(!x1) {
	printf("Memory allocation error");
		delete []A;
		delete []X;
		return -1;
    }
    x2 = new double[n];
    if(!x2) {
	printf("Memory allocation error");
		delete []A;
		delete []X;
		delete []x1;
		return -1;
    }

	for(int i = 0; i < n; i++) {
    	x1[i] = 0;
    	x2[i] = 0;
    }

    if(s != 0) {
        input_from_func(A, n, s);
    }
    else {
        name = argv[5];
        input = fopen(name, "r");
        if(!input) {
            printf("File open error");
            delete []A;
            delete []X;
            delete []x1;
            delete []x2;
            return -1;
        }
        fl_1 = input_from_file(A, n, input);
        if(fl_1 != 0) {
            printf("File reading error");
            fclose(input);
            delete []A;
            delete []X;
            delete []x1;
            delete []x2;
            return -1;
        }
        fclose(input);
    }


    printf("Matrix: \n");
    print_matrix(A, n, n, m);
    printf("\n");

    double norma = norm_matrix(n, A);

    pochti_triangle(n, A, x1); // Метод отражений

    //double *B = nullptr;
   // B = new double [n*n];

    //for (int s = 0; s < n; s++){
    //    for (int f = 0; f < n; f++)
    //        B[s*n+f] = A[f*n+s];}




    /*for (int s = 0; s < n; s++){
        for (int f = 0; f < n; f++)
            A[s*n+f] = B[f*n+s];}*/

    t2 = clock();
    int res = find_sobsn_znach(n, A, X, eps * norma, x1, x2); // Метод отражений
    t2 = (clock() - t2) / CLOCKS_PER_SEC;

        if(s != 0) {
            input_from_func(A, n, s);
        }
        else {
            input = fopen(name, "r");
            fl_1 = input_from_file(A, n, input);
            fclose(input);
        }
        if (!res)
          {

            r1 = Nevyazka_1(n, A, X);
            r2 = Nevyazka_2(n, A, X);
            printf("Eigenvalues: \n");
            print_matrix(X, 1, n, m);
            printf("\n");
          }
        else
          {
            r1 = -1.0;
            r2 = -1.0;

            printf ("No real Eigenvalues \n");
          }

    printf("Nevyazka_1 = %e Nevyazka_2 = %e Elapsed = %.2f secs\n", r1, r2, t2);

    delete []A;
    delete []X;
    delete []x1;
    delete []x2;
    return 0;
}
