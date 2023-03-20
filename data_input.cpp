#include <iostream>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "functions.h"
#include <climits>

double f(int k, int n, int i, int j)
{
  int row = i + 1;
  int col = j + 1;
  switch (k)
  {
  case 1:
    return n - std::max(row, col) + 1;
  case 2:
    return std::max(row, col);
  case 3:
    return fabs(row - col);
  case 4:
    return 1.0 / (row + col - 1);

  }

  return 0;
}


void print_matrix(int m, int l, int n, double *array)
{
  int row = m;
  int col = m;

  if (m > l)
    row = l;
  if (m > n)
    col = n;

  for (int i = 0; i < row; i++)
  {
      for (int j = 0; j < col; j++)
        printf(" %10.3e", array[i * l + j]);
    printf("\n");
  }
}


int data_input(double *A, double *b, int k, int n, FILE *input, int argc)
{
  double tmp = 0;

  if (argc == 4)
  {
    for (int i = 0; i < n; i++)
    {
      tmp = 0;
      for (int j = 0; j < n; j++)
      {
        A[i * n + j] = f(k, n, i, j);
        if (j % 2 == 0)
          tmp += A[i * n + j];
      }
      b[i] = tmp;
    }
  }
  else
  {

    for (int i = 0; i < n; i++)
    {
      tmp = 0;
      for (int j = 0; j < n; j++){
        if (fscanf(input, "%lf", &A[i * n + j]) != 1)
        {
          if (!feof(input))
          {
            fprintf(stdout, "Problems at element [%d,%d]\n",
                    i, j);
          }
          else
            fprintf(stdout, "Please, input more elements\n");
        return -1;
        }
        if (j % 2 == 0)
          tmp += A[i * n + j];
      }
      b[i] = tmp;
    }
  }

  return 1;
}


int correct_input(int argc, char **argv, char **file_name, int &n, int &m,
               int &k)
{
  if (argc == 4 || argc == 5)
  {
    if (argc == 4)
    {
      if (!((sscanf(argv[1], "%d", &n) == 1) && (n > 0) && (n < INT_MAX) && (sscanf(argv[2], "%d", &m) == 1) && (m >= 0) && (m < INT_MAX) && (sscanf(argv[3], "%d", &k) == 1) && (k >= 1 && k <= 4)))
      {

        return  -1;
      }
    }
    else
    {
      if (argc == 5)
      {
        if (!((sscanf(argv[1], "%d", &n) == 1) && (n > 0) && (n < INT_MAX) && (sscanf(argv[2], "%d", &m) == 1) && (m >= 0) && (m < INT_MAX) && (sscanf(argv[3], "%d", &k) == 1) && (k == 0)))
        {

          return  -1;
        }
        *file_name = argv[4];
        if (!file_name)
          return  -1;
      }
    }
  }
  else
  {

    return  -1;
  }

  return 1;
}




