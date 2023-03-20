#include <stdio.h>
#include <algorithm>
#include "functions.h"
#include <math.h>
#include "stdlib.h"


double f2(int i, int j) {
	if(i == j) {
		return 2;
	}
	else if(abs(i - j) == 1) {
		return -1;
	}
	else {
		return 0;
	}
}


double f3(int n, int i, int j) {
	if(i == j && j < n) {
		return 1;
	}
	else if(j == n) {
		return i;
	}
	else if(i == n) {
		return j;
	}
	else {
		return 0;
	}
}

double func(int s, int n, int i, int j) {
	int i_real = i + 1;
	int j_real = j + 1;
    switch(s) {
        case 1:
            return n - std::max(i_real, j_real) + 1;
        case 2:
            return f2(i_real, j_real);
        case 3:
            return f3(n, i_real, j_real);
        case 4:
            return 1/((double)(i_real + j_real - 1));
    }
    return -1;
}


void input_from_func(double *A, int n, int s) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            A[i * n + j] = func(s, n, i, j);
        }
    }
}


int input_from_file(double *A, int n, FILE *input) {
    double elem = 0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(fscanf(input, "%lf", &elem) != 1){
                printf("input_from_file: failed to read the next element\n");
                return -1;
            }
            A[i * n + j] = elem;
        }
    }
    return 0;
}


void print_matrix(double *A, int l, int n, int r) {
    int rows = 0, cols = 0;
    rows = (r < l? r : l);
    cols = (r < n? r : n);
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            printf("%10.3e ", A[i * n + j]);
        }
        printf("\n");
    }
}

