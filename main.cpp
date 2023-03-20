#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include "time.h"
#include <math.h>


int main(int argc, char **argv)
{
  FILE *f  = 0;
  char *file_name;

  int n = 0, r = 0, k = 0, fl_1 = 0, fl_2 = 0, fl_solve = 0;
  double seconds;

  double *A, *b, *x;
  int *index;

  clock_t t1, t2;

  fl_1 = correct_input(argc, argv, &file_name, n, r, k);

  if (fl_1 == -1){
    std::cout << "Data initialization error";
    return -1;}

  if (argc == 5)
  {
    f = fopen(file_name, "rt");
    if (f == 0)
    {
      std::cout << "File open error";
      return -1;
    }
  }

  A = new double[n * n];
  if (A == NULL)
  {
    std::cout << "Memory allocation error";
    if (argc == 5)
    fclose(f);
    return -1;
  }
  b = new double[n];
  if (b == NULL)
  {
    delete[] A;
    std::cout << "Memory allocation error";
    if (argc == 5)
    fclose(f);
    return -1;
  }
  fl_2 = data_input(A, b, k, n, f, argc);

  if (fl_2 == -1)
  {
    delete[] A;
    delete[] b;
    if (argc == 5)
      fclose(f);
    std::cout << "Data initialization error";
    return -1;
  }


  std::cout<<"\nMatrix:\n";
  print_matrix(r, n, n, A);


  x = new double[n];
  if (x == NULL)
  {
    delete[] A;
    delete[] b;
    if (argc == 5)
    fclose(f);
    return -1;
  }

  index = new int[n];
  if (index == NULL)
  {
    delete[] A;
    delete[] b;
    delete[] x;
    if (argc == 5)
    fclose(f);
    return -1;
  }
  for (int i = 0; i < n; i++)
    index[i] = i;

  t1 = clock();
  fl_solve = Jordan_mid(n, A, b, x, index);
  t2 = clock ();
  seconds = double(t2 - t1)/CLOCKS_PER_SEC;

  if (fl_solve == -1){
    std::cout<<"\nInvalid matrix";
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] index;
    if (argc == 5)
    fclose(f);
    return -1;}

  std::cout<<"\nSolution:\n";
  print_matrix(r, 1, n, x);

  printf("\nThe time: %.2f seconds\n", seconds );

  std::cout<<"Sin_norm: ";

  std::cout<<sin_norm(x, n)<<std::endl;

  if (argc == 5)
    rewind(f);

  fl_2 = data_input(A, b, k, n, f, argc);
  std::cout<<"Nevyazka:";
  printf("%10.3e\n", Nevyazka(A, x, b, n));

  if (argc == 5)
    fclose(f);

  delete[] A;
  delete[] b;
  delete[] index;
  delete[] x;
  return 0;
}
