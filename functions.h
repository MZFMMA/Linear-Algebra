
double func(int s, int n, int i, int j);

void input_from_func(double *A, int n, int s);

int input_from_file(double *A, int n, FILE *input);

void print_matrix(double *A, int n, int l, int r);

void pochti_triangle(int n, double *A, double *x);

double Nevyazka_1(int n, double *A, double *X);

double Nevyazka_2(int n, double *A, double *X);

double norm_matrix(int n, double *A);

int find_sobsn_znach(int n, double *A, double *X, double eps, double *x1, double *x2);



void almost_triangle(int n, double *A, double *x);
