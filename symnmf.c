#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "symnmf.h"


double** create_zero_mat(int n, int d)
{/* return a zero matrix order n x d */
    double **mat; 
    int i;
    mat = (double**) calloc(n, sizeof(double*));
    if (NULL == mat)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        mat[i] = (double*) calloc(d, sizeof(double));
        if (NULL == mat[i])
        {
            return NULL;
        }
    }
    return mat;
}

void freeMat(double **mat, int n)
{ /*free memory of matrice*/
    int i;
    for(i = 0; i < n; i++){
        free(mat[i]);
    }
    free(mat);
}

double** arr_cast_mat(double *arr, int n, int d)
{ /*Take an array and transorm it into matrice*/
    int i,j;
    double **mat;
    mat = create_zero_mat(n,d);
    if (mat == NULL)
    {
        printf("An Error Has Occurred");
        return NULL;
    }
    for(i = 0; i < n ;i++)
    {
        for(j = 0; j < d; j++){
            mat[i][j] = arr[i * d + j];
        }
    }
    return mat;
}

double** create_A_mat(double **points, int n, int d)
{/* return A - Simlarity matrice  */
    double **A_mat;
    double delta = 0;
    double value;
    int i,j, k;
    double sum = 0;

    A_mat = create_zero_mat(n,n);

    if (NULL == A_mat)
    {
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if (i == j){
                A_mat[i][j] = 0;
            }
            else
            {
                sum = 0;
                for(k = 0; k < d ;k++)
                {
                    delta = pow((points[i][k] - points[j][k]), 2);
                    sum += delta;
                }
                value = ((-1) * sum) / 2;
                A_mat[i][j] = exp(value);
            }
        }
    }
    return A_mat;
}

double** create_D_mat(double **A_mat, int n)
{/* return D - Diagonal Degree Matrix */
    double **D_mat;
    double sum;
    int i ,j;

    D_mat = create_zero_mat(n, n);

    if (NULL == D_mat)
    {
        return NULL;
    }

    sum = 0.0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            sum = sum + A_mat[i][j];
        }
        D_mat[i][i] = sum;
        sum = 0.0;
    }
    return D_mat;
}

double** D_mat_normalized(double **A_mat, int n)
{/* return D -  Diagonal Degree Matrix ^(-1/2) */
    double **D_mat;
    double sum;
    int i ,j;

    D_mat = create_zero_mat(n, n);

    if (NULL == D_mat)
    {
        return NULL;
    }

    sum = 0.0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            sum = sum + A_mat[i][j];
        }
        D_mat[i][i] = 1/sqrt(sum);
        sum = 0.0;
    }
    return D_mat;
}

double** mat_multipication(double **matA, double **matB, int n)
{   /* return an array that is the result of A*B matrices */ 
    double **res;
    int j;
    int k;
    int i;
    double sum;

    res = create_zero_mat(n,n);
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++)
        {
            sum = 0;
            for (k = 0; k < n; k++)
            {
                sum += matA[i][k]*matB[k][j];
            }
            res[i][j]= sum;
        }
    }
    return res;
}

double** create_W_mat(double **points, int n, int d)
{
    /* return W - the normalized similarity matrix - W = D^(-1/2) * A * D^(-1/2) */
    double **A_mat, **D_mat;
    double **first_multipication, **second_multipication;

    A_mat = create_A_mat(points, n, d);
    if (NULL == A_mat)
    {
        return NULL;
    }
    
    D_mat = D_mat_normalized(A_mat, n);
    if (NULL == D_mat )
    {
        return NULL;
    }

    first_multipication = mat_multipication(D_mat, A_mat, n);
    if (NULL == first_multipication)
    {
        return NULL;
    }
    second_multipication = mat_multipication(first_multipication, D_mat, n);
    if (NULL == second_multipication )
    {
        return NULL;
    }

    freeMat(A_mat, n);
    freeMat(D_mat, n);
    freeMat(first_multipication, n);
    return second_multipication;
}


void printm(double **mat, int n, int d)
{
    /* prints matrice  */
    int i,j;
    for ( i = 0; i < n; i++)
    {
        for ( j = 0; j < d; j++)
        {
            printf("%.4f", mat[i][j]);
            if (j < (d - 1))
            {
                printf(",");
            }
            else{
                printf("\n");
            }
        }
    }
}

double** transpose(double **mat, int n){
    /* return transpose of matrice */
    double **mat_T;
    int i,j;
    mat_T = create_zero_mat(n, n);
    if (mat_T == NULL){
        return NULL;
    }
    for(i = 0; i < n; i++){
        for( j = 0; j < n; j++){
            mat_T[j][i] = mat[i][j];
        }
    }
    return mat_T;
}

double frobenius_norm(double **mat, int n, int k){
    /* return frobenius norm of matrice */
    int i,j;
    double norm = 0;
    for(i = 0; i < n; i++){
        for( j = 0; j < k; j++){
            norm += pow(mat[i][j], 2);
        }
    }
    return norm;
}

void get_points_input(FILE *ifp, double* points_arr)
{/* points from the file to an array */
    double point = 0.0;
    int i=0;
    while (fscanf_s(ifp, "%lf,", &point) != EOF)
    {
        points_arr[i] = point;
        i++;
    }
}
double** update_H_mat(double **matH, double **matW, int n, int k, double epsilon, int iter){
    double **new_H, **numerator_H, **denominator_H, **mat_H_T, **frobenius_mat;
    int i,j;
    double norm;
    
    numerator_H = mat_multipication(matW, matH, n);
    mat_H_T = transpose(matH, n);
    denominator_H = mat_multipication(matH, mat_multipication(mat_H_T, matH, n), n);
    new_H = create_zero_mat(n, k);
    frobenius_mat = create_zero_mat(n, k);

    for (i = 0; i < n; i++){
        for (j = 0; j < k; j++){
            new_H[i][j] = matH[i][j] * (0.5 + ((0.5 * numerator_H[i][j]) / denominator_H[i][j]));
            frobenius_mat[i][j] = new_H[i][j] - matH[i][j];
        }
    }

    norm = frobenius_norm(frobenius_mat, n, k);
    if (iter == 0 || norm < epsilon){
        return matH;
    }
    return update_H_mat(new_H, matW, n, k, epsilon, iter);
}

int main(int argc, char** argv){
    /* if there are command-line arguments, they are interpered as filenames, and processed in order */
    FILE *ifp; 
    int d; 
    int n;  
    char *filename;
    double *data_points;
    double **points;
    double **A_mat, **D_mat, **W_mat;
    char charCount;
    int err;

    filename = argv[2];
    err = fopen_s(&ifp, filename, "r");
    if (err != 0){
        printf("An Error Has Occurred");
        return 1;
    }

    /* check if function ok */
    n=0;
    d = 1;
    while ((charCount = getchar()) != EOF) /* switch with fget? */
    {
        if (charCount == '\n') {
            n = n+1;
        } 
        else{
            if (n == 0 && charCount == ',') 
            {
            d = d+1;
            } 
        }
    }
    rewind(stdin); /* reset pointer*/

    data_points = calloc(n*d, sizeof(double));
    if (NULL == data_points){
        printf("An Error Has Occurred");
        return 1;
    }
    get_points_input(ifp,data_points);
    fclose(ifp);

    points = arr_cast_mat(data_points,n,d);
    if (points == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    if (strcmp(argv[1],"sym") == 0)
    {
        A_mat = create_A_mat(points, n, d);
        if (NULL == A_mat)
        {
            printf("An Error Has Occurred");
            return 1;
        }
        printm(A_mat, n, n);
        freeMat(A_mat, n);
        freeMat(points,n);
    }

    if (strcmp(argv[1],"ddg") == 0)
    {
        D_mat = create_D_mat(points, n);
        if (NULL == D_mat)
        {
            printf("An Error Has Occurred");
            return 1;
        }
        printm(D_mat, n, n);
        freeMat(D_mat, n);
        freeMat(points,n);
    }

    if (strcmp(argv[1],"norm") == 0)
    {
        W_mat = create_W_mat(points, n, d);
        if (NULL == W_mat)
        {
            printf("An Error Has Occurred");
            return 1;
        }
        printm(W_mat, n, n);
        freeMat(W_mat, n);
        freeMat(points,n);
    }

    if (strcmp(argv[1],"H") == 0)
    {
        W_mat = create_W_mat(points, n, d);
        if (NULL == W_mat)
        {
            printf("An Error Has Occurred");
            return 1;
        }
        printm(W_mat, n, n);
        freeMat(W_mat, n);
        freeMat(points,n);
    }
    free(data_points);
    return 0;
}
