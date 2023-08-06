# ifndef symnmfh
# define symnmfh

double** create_zero_mat(int n, int d);
void freeMat(double **mat, int n);
double** arr_cast_mat(double *arr, int n, int d);
double** create_A_mat(double **points, int n, int d);
double** create_D_mat(double **A_mat, int n);
double** D_mat_normalized(double **A_mat, int n);
double** mat_multipication(double **matA, double **matB, int n);
double** create_W_mat(double **points, int n, int d);
void printm(double **mat, int n, int d);
double** transpose(double **mat, int n);
double frobenius_norm(double **mat, int n, int k);
void get_points_input(FILE *ifp, double* points_arr);
double** update_H_mat(double **matH, double **matW, int n, int k, double epsilon, int iter);


# endif