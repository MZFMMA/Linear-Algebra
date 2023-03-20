double matrix_norm(double *matrix, int n);
int Jordan_mid(int n, double* a, double* b, double* x, int* index);
int data_input(double *A, double *b, int k, int n, FILE *input, int argc);
double f(int k, int n, int i, int j);
int correct_input(int argc, char **argv, char **file_name, int &n, int &r, int &k);
void print_matrix(int m, int l, int n, double *array);
double Nevyazka(double *A, double *x, double *b, int n);
double sin_norm(double *x, int n);