# ifndef symnmfh
# define symnmfh
# define _CRT_SECURE_NO_WARNINGS 1

double** create_zero_mat(int n, int d);
void freeMat(double **mat, int n);
double** arr_cast_mat(double *arr, int n, int d);
double** create_A_mat(double **points, int n, int d);
double** create_D_mat(double **A_mat, int n);
double** D_mat_normalized(double **A_mat, int n);
void mat_multipication(double **matA, double **matB, double **matC, int n, int m, int d);
double** create_W_mat(double **points, int n, int d);
void printm(double **mat, int n, int d);
void transpose(double **mat, double **matT, int n, int d);
double frobenius_norm(double **matA,double **matB, int n, int k);
double** get_points_input(FILE *ifp, int n, int d);
double** update_H_mat(double **matH, double **matW, int n, int k, double epsilon, int iter);
double **calcH(double** matH, double** matW, int n, int k);

# endif