#include <stdio.h>
#include <math.h>
#include <stdlib.h>



double** create_zero_mat(int n, int d)
{
    /* return a zero matrix order n x d */
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
{
    int i;
    for(i = 0; i < n; i++){
        free(mat[i]);
    }
    free(mat);
}

double** create_A_mat(double **points, int n, int d)
{
    /* return A - Simlarity matrice  */
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
{
    /* return D - Diagonal Degree Matrix */
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
{
    /* return D -  Diagonal Degree Matrix ^(-1/2) */
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
{
    
    double **res;
    int j;
    int k;
    int i;
    double sum;

    res = zero_mat(n,n);
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

double** W_mat_laplacian(double **points, int n, int d)
{
    /* return W - the normalized similarity matrix - D^(-1/2)AD^(-1/2) */
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

    first_multipication = mat_mult(D_mat, A_mat, n);
    if (NULL == first_multipication)
    {
        return NULL;
    }
    second_multipication = mat_mult(first_multipication, D_mat, n);
    if (NULL == second_multipication )
    {
        return NULL;
    }

    freeMat(A_mat, n);
    freeMat(D_mat, n);
    freeMat(first_multipication, n);
    return second_multipication;
}

double** transpose(double **mat, int n)
{
    /* return mat^t - the transpose of the input matrice */
    double **mat_T;
    double temp;
    int i,j;
    mat_T = create_zero_mat(n, n);
    if (NULL == mat_T)
    {
        return NULL;
    }
    for(i = 0; i < n; i++)
    {
        for( j = 0; j < n; j++)
        {
            temp = mat[i][j];
            mat_T[j][i] = temp;
        }
    }
    return mat_T;
}

double avg_mat_entries(double **W_mat, int n)
{
    /* return m - the average value of all entries of W */
    int i,j;
    double m;
    if(NULL == n)
    {
        reutrn NULL;
    }
    for(i = 0; i < n; i++)
    {
        for( j = 0; j < n; j++)
        {
            m += W_mat[i][j];

        }
    }
    m = m / (n * n);

    return m;
}

double** create_H_Mat(double **W_mat, int n ,int k)
{   
    double **H_mat;
    int i,j,endpoint,temp;
    double m; 
    m=avg_mat_entries( W_mat , n );
    endpoint = 2*sqrt(m\k);
    
    H_mat =create_zero_mat(n,k);
    for(i = 0; i < n; i++)
    {
        for( j = 0; j < n; j++)
        {
            temp=( rand() ) % (endpoint+1); /* get number in range [0,endpoint]*/
            H_mat[i][j] = temp;

        }
    }

    return H_mat;
}