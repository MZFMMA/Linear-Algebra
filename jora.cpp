#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include "functions.h"

double matrix_norm(double *matrix, int n)
{
  double norm = 0, tmp_sum = 0;
  for (int i = 0; i < n; i++)
  {
    tmp_sum = 0;
    for (int j = 0; j < n; j++)
      tmp_sum += fabs(matrix[i * n + j]);
    if (tmp_sum > norm)
      norm = tmp_sum;
  }
  return norm;
}

double Nevyazka(double *A, double *x, double *b, int n)
{
    double *Axb;
    Axb = new double [n];
    for (int i=0; i<n; i++)
        Axb[i]=0;
    double n1 = 0, n2 = 0;

    for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    {
                    Axb[i] += A[i * n + j] * x[j];

                    }
                Axb[i] -= b[i];
                n1 += Axb[i]*Axb[i];

                n2 += b[i]*b[i];
            }


    n1 = sqrt(n1);
    n2 = sqrt(n2);
    delete[] Axb;
    return n1/n2;
}

double sin_norm(double *x, int n)
{
    double sum = 0;
    int t;
    for(int i = 0; i < n; i++)
    {
        if (i % 2 == 0)
            t = 1;
        else
            t = 0;
        sum += pow (x[i] - t, 2);
    }

    return sqrt(sum);
}





int Jordan_mid(int n, double* A, double* b, double* x, int* index)
{
	int max_row, max_col, tr;
	double tmp, max_el;


    double norm = matrix_norm(A, n);

	for (int i = 0; i < n; i++)
		index[i] = i;

	for (int i = 0; i < n; i++)
	{
		max_el = fabs(A[i * n + i]);
		max_row = i;
		max_col = i;

		for ( int j = i; j < n; j++)
			for (int k = i; k < n; k++)
				if (max_el < fabs(A[j * n + k]))
				{
					max_el = fabs(A[j * n + k]);
					max_row = j;
					max_col = k;
				}


                if (fabs(A[max_row * n + max_col]) <= 1.2e-16 * norm)

                    return -1;

		tr = index[i];
		index[i] = index[max_col];
		index[max_col] = tr;

		tmp = b[i];
		b[i] = b[max_row];
		b[max_row] = tmp;

		for (int j = 0; j < n; j++)
		{
			tmp = A[i * n + j];
			A[i * n + j] = A[max_row * n + j];
			A[max_row * n + j] = tmp;
		}

        for (int j = 0; j < n; j++)
		{
			tmp = A[j * n + i];
			A[j * n + i] = A[j * n + max_col];
			A[j * n + max_col] = tmp;
		}


        tmp = 1.0 / A[i * n + i];
		for (int j = i; j < n; j++)
			A[i * n + j] *= tmp;
		b[i] *= tmp;

		for (int j = 0; j < n; j++)
		{
			if (j<i)
			{tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}

            if (j>i){
			tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}

			else{}
		}

                //std::cout<<"\n--------------------\n";
                //print_matrix(n,n,n, A);

        }



	for (int i = 0; i < n; i++)
		x[index[i]] = b[i];

	return 0;}

/*int main(void)
{
double a[]= {3, 2, -5, 2, -1, 3, 1, 2, -1};
double b[] = {-1, 13 , 9};
double x[] = {0, 0, 0};
int index[] = {0, 0, 0, 0, 0, 0};
int ara;
ara = Jordan_mid(3, a, b, x ,index);
std::cout<<ara<<"\n";
std::cout<<x[0]<<"\n";
std::cout<<x[1]<<"\n";
std::cout<<x[2]<<"\n";
return 1;
}*/


