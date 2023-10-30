# ifndef symnmfh
# define symnmfh
# define _CRT_SECURE_NO_WARNINGS 1

double** create_zero_mat(int n, int d);
void freeMat(double **mat, int n);
double** create_A_mat(double **points, int n, int d);
double** create_D_mat(double **points, int n, int d);
double** D_mat_normalized(double **points, int n, int d);
void mat_multipication(double **matA, double **matB, double **matC, int n, int m, int d);
double** create_W_mat(double **points, int n, int d);
void printm(double **mat, int n, int d);
void transpose(double **mat, double **matT, int n, int d);
double frobenius_norm(double **matA,double **matB, int n, int k);
double** update_H_mat(double **matH, double **matW, int n, int k, double epsilon, int iter);
double **calcH(double** matH, double** matW, int n, int k);
int* getDimensions(char* file_name);
double** getPoints(char* file_name, int n, int d);

# endif
