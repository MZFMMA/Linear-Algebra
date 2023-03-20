#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include "functions.h"
#include <sys/time.h>


void
reduce_sum(int p, double *a = nullptr, int n = 0);
void reduce_sum(int p, double *a, int n)
{
  static pthread_mutex_t m     = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t  c_in  = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t  c_out = PTHREAD_COND_INITIALIZER;
  static int t_in  = 0;
  static int t_out = 0;
  static double *r = nullptr;
  int i;

  if (p <= 1)
    return;
  pthread_mutex_lock(&m);
  if (r == nullptr)
    r = a;
  else
    for (i = 0; i < n; i++)
      r[i] += a[i];
  t_in++;
  if (t_in >= p)
  {
    t_out = 0;
    pthread_cond_broadcast(&c_in); // оповещение что все потоки вошли
  }
  else
    while (t_in < p)
      pthread_cond_wait(&c_in, &m);
  if (r != a)
    for (i = 0; i < n; i++)
      a[i] = r[i];
  t_out++;
  if (t_out >= p)
  {
    t_in = 0;
    r = nullptr;
    pthread_cond_broadcast(&c_out);
  }
  else
    while (t_out < p)
      pthread_cond_wait(&c_out, &m);
  pthread_mutex_unlock(&m);
}


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
                n1 += pow(Axb[i],2);

                n2 += pow(b[i],2);
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

//-------------------------------------------------------------------------

struct thread_data{
    double* A;
    double* b;
    double* t;
	int* indicies;
	int n;
	int index;
	sem_t* sem_1;
	sem_t* sem_2;
	pthread_mutex_t* mutex;
	int* fl;
	int* thread_finished_counter;
	int number_of_threads;
};

int Jora_parallel(int n, double* A, double* b, double* x, int* indicies, int number_of_threads, double *t){

    pthread_t *threads;
    struct thread_data *threads_data;


    int fl_2 = 0;


    sem_t sem_1;
	sem_init(&sem_1, 0, 0);

	sem_t sem_2;
	sem_init(&sem_2, 0, 0);

    int fl = 1;
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex, NULL);

       threads = new pthread_t[number_of_threads];
       // pthread_t threads[number_of_threads];
	int thread_finished_counter = 0;

        threads_data = new struct thread_data[number_of_threads];
      //  struct thread_data threads_data[number_of_threads];
	for (int i = 0; i < number_of_threads; i++){
		threads_data[i].A = A;
		threads_data[i].b = b;
        threads_data[i].t = t;
		threads_data[i].n = n;
        threads_data[i].indicies = indicies;
		threads_data[i].sem_1 = &sem_1;
		threads_data[i].fl = &fl;
		threads_data[i].sem_2 = &sem_2;
		threads_data[i].mutex = &mutex;
		threads_data[i].thread_finished_counter = &thread_finished_counter;
		threads_data[i].index = i;
		threads_data[i].number_of_threads = number_of_threads;
	}

	for (int i = 1; i < number_of_threads; i++){
                if((pthread_create(&threads[i], NULL, Jora_thread_func,(void*) &threads_data[i]))!=0)
                {
                    printf("%d",i);
                    fl_2 = -1;}
                //std::cout<<"\nN - "<< thread_finished_counter<<std::endl;

	}

        if (fl_2 == -1)
        {
            delete [] threads;
            delete [] threads_data;
            return -2;
        }


    Jora_thread_func(threads_data + 0);
	for (int i = 1; i < number_of_threads; i++){

		if((pthread_join(threads[i], NULL))!=0)
        {
            delete [] threads;
            delete [] threads_data;
              printf("efwf");
            return -1;
        }
        //std::cout<<"\nM - "<< thread_finished_counter<<std::endl;
	}


    if (fl == -69){
        delete [] threads;
        delete [] threads_data;
        return -1;}

    for (int i = 0; i < n; i++)
		x[indicies[i]] = b[i];
    delete [] threads;
    delete [] threads_data;

return 0;
}

void* Jora_thread_func(void* arg){

    struct thread_data* data = (struct thread_data*)arg;
    int max_row, max_col, tr;
    double tmp, max_el;
    double* A = data->A;
    double* b = data->b;

	int* indicies = data->indicies;
	int n = data->n;
   	int p = data->index;
    int N = data->number_of_threads;
    double norm = matrix_norm(A, n);


    struct timeval t1, t2;

    clock_t tim;
    reduce_sum(N);
    gettimeofday(&t1, NULL);

    tim = clock();

	for (int i = 0; i < n; i++) {


		if(p == 0)
        {

            max_el = fabs(A[i * n + i]);
            max_row = i;
            max_col = i;

            for (int j = i; j < n; j++)
                for (int k = i; k < n; k++)
                    if (max_el < fabs(A[j * n + k]))
                    {
                        max_el = fabs(A[j * n + k]);
                        max_row = j;
                        max_col = k;
                    }

            if (fabs(A[max_row * n + max_col]) <= 1.2e-16 * norm){

                *data->fl = -69;


                return NULL;}

            tr = indicies[i];
            indicies[i] = indicies[max_col];
            indicies[max_col] = tr;

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

        }


        reduce_sum(N);
        if (*data->fl == - 69)
            return NULL;





		for (int j = 0 + p; j < n; j+=N)
        {

			if (j!=i)
			{tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}

            /*if (j>i){
			tmp = A[j * n + i];
			for (int k = i; k < n; k++)
				A[j * n + k] -= A[i * n + k] * tmp;
			b[j] -= b[i] * tmp;}*/


        }

        reduce_sum(N);
        }
    //std::cout << p <<"\tTIME WORKING:\t"<< clock() - tim<< std::endl;

    gettimeofday(&t2, NULL);
    *data->t = (t2.tv_sec - t1.tv_sec)*1000;
    *data->t += (t2.tv_usec - t1.tv_usec)/1000;
    *data->t /= 1000;


	return NULL;
}






