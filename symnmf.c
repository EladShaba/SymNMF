#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "symnmf.h"

/* struct for dots' matrix size (n x d)*/
typedef struct dot_size{
   int n;
   int d;
}SIZE;

/* create a n x d zero matrice */
double** create_zero_mat(int n, int d){
    double **mat; 
    int i;
    
    mat = (double**) calloc(n, sizeof(double*));
    if (mat == NULL){
        printf("An Error Has Occurred");
        return NULL;
    }

    for (i = 0; i < n; i++){
        mat[i] = (double*) calloc(d, sizeof(double));
        if (mat[i] == NULL){
            printf("An Error Has Occurred");
            return NULL;
        }
    }

    return mat;
}

/* free memory of matrice */
void freeMat(double **mat, int n){ 
    int i;
    for(i = 0; i < n; i++)
        free(mat[i]);
    free(mat);
}

/* create and return the A matrix - Simlarity matrice */
double** create_A_mat(double **points, int n, int d){    
    double **A_mat, value, sum;
    int i, j, k;

    if (points == NULL)
        return NULL;

    A_mat = create_zero_mat(n,n);
    if (A_mat == NULL)
        return NULL;

    for (i = 0; i < n; i++)
        for(j = 0; j < n; j++){
            if (i != j){
                sum = 0;
                for(k = 0; k < d ;k++)
                    sum += pow((points[i][k] - points[j][k]), 2);
                value = ((-1)*sum) / 2.0;
                A_mat[i][j] = exp(value);
            }
        }
    
    return A_mat;
}

/* create and return the D matrix - Diagonal Degree Matrix */
double** create_D_mat(double **A_mat, int n){
    double **D_mat, sum;
    int i ,j;

    if (A_mat == NULL){
        return NULL;
    }
    D_mat = create_zero_mat(n, n);
    if (D_mat == NULL){

    
        return NULL;
    }
    for (i = 0; i < n; i++){
        sum = 0.0;
        for (j = 0; j < n; j++)
            sum += A_mat[i][j];
        D_mat[i][i] = sum;
    }
    
    return D_mat;
}

/* create and return the normalized D matrix - Diagonal Degree Matrix (D^(-1/2)) */
double** D_mat_normalized(double **A_mat, int n){
    double **D_mat;
    int i;
    
    if (A_mat == NULL)
        return NULL;
    
    D_mat = create_D_mat(A_mat, n);
    if (D_mat == NULL)
        return NULL;

    for (i = 0; i < n; i++)
        D_mat[i][i] = 1 / sqrt(D_mat[i][i]);

    return D_mat;
}

/* fill C with the result of A*B where
A is a n x m matrice, B is a m x d matrice and C is a n x d matrice */ 
void mat_multipication(double **matA, double **matB, double **matC, int n, int m, int d){   
    double sum;
    int i, j, k;

    for (i = 0; i < n; i++)
        for (j = 0; j < d; j++){
            sum = 0.0;
            for (k = 0; k < m; k++)
                sum += matA[i][k] * matB[k][j];
            matC[i][j]= sum;
        }
}

/* create and return the W matrix - the normalized similarity matrix - W = D^(-1/2) * A * D^(-1/2) */
double** create_W_mat(double **points, int n, int d){
    double **A_mat, **D_mat;
    double **first_multipication, **second_multipication;
    
    first_multipication = create_zero_mat(n, n);
    second_multipication = create_zero_mat(n, n);
    if (first_multipication == NULL || second_multipication == NULL)
        return NULL;
    
    A_mat = create_A_mat(points, n, d);
    D_mat = D_mat_normalized(A_mat, n);
    if (A_mat == NULL || D_mat == NULL)
        return NULL;
    mat_multipication(D_mat, A_mat, first_multipication, n, n, n);
    mat_multipication(first_multipication, D_mat, second_multipication, n, n, n);

    /*free memory*/
    if (A_mat != NULL)
        freeMat(A_mat, n);
    if (D_mat != NULL)
        freeMat(D_mat, n);
    if (first_multipication != NULL)
        freeMat(first_multipication, n);

    return second_multipication;
}

/* prints a n x d matrice */
void printm(double **mat, int n, int d){
    int i,j;
    
    for ( i = 0; i < n; i++){
        for ( j = 0; j < d - 1; j++)
            printf("%.4f,", mat[i][j]);
        printf("%.4f\n", mat[i][d-1]);
    }
}

/* parse matT such that matT is the transposed matrice of mat (which is a n x d matrice) */
void transpose(double **mat, double **matT, int n, int d){
    int i,j;
    
    for(i = 0; i < n; i++)
        for( j = 0; j < d; j++)
            matT[j][i] = mat[i][j];
}

/* return frobenius norm of matrice */
double frobenius_norm(double **matA,double **matB, int n, int k){
   double dist =0.0;
    int i,j;
    for(i=0;i<n;i++)
    {

    for(j=0; j<k; j++)
        {
        dist += (matA[i][j]-matB[i][j])*(matA[i][j]-matB[i][j]);
        }
    }
    return dist;

}



double **calcH(double** matH, double** matW, int n, int k){
    int i,j;
    double B = 0.5;
    double **newH,** matWH,**matHH,**matHHH,**matHt;
    newH=create_zero_mat(n,k);
    matWH=create_zero_mat(n,k);
    matHH=create_zero_mat(n,n);
    matHHH=create_zero_mat(n,k);
    matHt=create_zero_mat(k,n);
    
    mat_multipication(matW,matH,matWH,n,n,k);
    transpose(matH,matHt,n,k);
    mat_multipication(matH,matHt,matHH,n,k,n);
    mat_multipication(matHH,matH,matHHH,n,n,k);
   
    for(i=0;i<n;i++)
    {
        for(j=0;j<k;j++){
            newH[i][j] = matH[i][j]*(1 -B + B*(matWH[i][j]/matHHH[i][j]));
        }
    }

    freeMat(matWH,n);
    freeMat(matHH,k);
    freeMat(matHHH,n);
    freeMat(matHt,k);
    return newH;
}

/* return the updated H matrix */
double** update_H_mat(double **matH, double **matW, int n, int k, double eps, int iter){
    int count = 0;
    double **newH;
    double **copyH;
    copyH =  matH;
    newH=calcH(copyH,matW,n,k);
    while( frobenius_norm(copyH,newH,n,k)>= eps && count < iter){
        copyH = newH;
        newH = calcH(copyH,matW,n,k);
        count++;
    }
    freeMat(copyH,n);
    return newH;
}
/* reading the file */
SIZE* reading_file(char* file_name, double*** dots){
    FILE *f = fopen(file_name,"r");
    int d = 0;
    int i,j;
    int row = 0;
    int col = 0;
    char c = 0;
    double num = 0;
    SIZE *ds = malloc(sizeof(SIZE));
    if (NULL == ds) { /* malloc check */
        return NULL;
    }

    /* reading the file the 1st time to calculate n and d*/
    while (fscanf(f,"%lf",&num)!= EOF){
        c = fgetc(f); /* reading the char after the number */
        if(EOF != c){
            /* same point */
            if(',' == c){
                col++;
            }
            else if ('\n' == c){
                row++;
                d = col;
                col = 0;
            }
        }
    }
    d++;
    ds->n = row;
    ds->d = d;
    
    /* reading the file the 2nd time enter dots to the array*/
    *dots = malloc(row * sizeof(double*));
    if (NULL == *dots) { /* malloc check */
        return NULL;
    }
    fseek(f, 0, SEEK_SET); /* return to the start of the file*/
    for(i = 0; i < row; i++){
        (*dots)[i] = malloc(d*sizeof(double)); 
        if (NULL == (*dots)[i]) { /* malloc check */
            return NULL;
        }
        for(j = 0; j < d; j++){
            if(fscanf(f,"%lf",&num)!= EOF){
                (*dots)[i][j] = num;
            }
            c = fgetc(f);
        }
    }

    fclose(f);
    return ds; /* return row = number of dots */
}
int main(int argc, char** argv){
    /* if there are command-line arguments, they are interpered as filenames, and processed in order */
    SIZE *ds=NULL;
    int n, d;   
    char *filename;
    
    double **points, **A_mat, **D_mat, **W_mat;
    if(argc != 3)
    {
        printf("An Error Has Occurred");
        return 1;
    }
   /* reading the file and calculate n and d */
    
    filename = argv[2];
    ds = reading_file(filename, &points);
    if (NULL == ds) {
        puts("An Error Has Occurred");
        return 1;
    }
    n = ds->n;
    d = ds->d;
    free(ds);

    if (strcmp(argv[1], "sym") == 0){
        A_mat = create_A_mat(points, n, d);
        if (A_mat == NULL)
            return 1;
        printm(A_mat, n, n);
        freeMat(A_mat, n);
    }

    if (strcmp(argv[1], "ddg") == 0){
        A_mat = create_A_mat(points, n, d);
        if (A_mat == NULL)
            return 1;
        D_mat = create_D_mat(A_mat, n);
        if (D_mat == NULL)
            return 1;
        
        printm(D_mat, n, n);
        freeMat(A_mat, n);
        freeMat(D_mat, n);
    }

    if (strcmp(argv[1], "norm") == 0){
        W_mat = create_W_mat(points, n, d);
        if (W_mat == NULL)
            return 1;
        
        printm(W_mat, n, n);
        freeMat(W_mat, n);
    }

    freeMat(points, n);
    return 0;
}