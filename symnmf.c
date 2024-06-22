#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "symnmf.h"


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

/* create and return Simlarity matrice */
double** create_A_mat(double **points, int n, int d){    
    double **A_mat, value, sum;
    int i, j, k;

    A_mat = create_zero_mat(n,n);
    if (A_mat == NULL)
        return NULL;

    for (i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if (i != j){
                sum = 0;
                for(k = 0; k < d ;k++)
                    sum += pow((points[i][k] - points[j][k]), 2);
                value = ((-1)*sum) / 2.0;
                A_mat[i][j] = exp(value);
            }
        }
    }

    return A_mat;
}

/* create and return Diagonal Degree Matrix */
double** create_D_mat(double **points, int n, int d){
    double **D_mat, sum, **A_mat;
    int i, j;
    
    A_mat = create_A_mat(points, n, d);
    if (A_mat == NULL)
        return NULL;

    D_mat = create_zero_mat(n, n);
    if (D_mat == NULL){
        freeMat(A_mat, n);
        return NULL;
    }

    for (i = 0; i < n; i++){
        sum = 0.0;
        for (j = 0; j < n; j++)
            sum += A_mat[i][j];
        D_mat[i][i] = sum;
    }

    freeMat(A_mat, n);
    return D_mat;
}

/* create and return Diagonal Degree Matrix (D^(-1/2)) */
double** D_mat_normalized(double **points, int n, int d){
    double **D_mat;
    int i;

    D_mat = create_D_mat(points, n,d);
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

    for (i = 0; i < n; i++){
        for (j = 0; j < d; j++){
            sum = 0.0;
            for (k = 0; k < m; k++)
                sum += matA[i][k] * matB[k][j];
            matC[i][j]= sum;
        }
    }
}

/* create and return the W matrix - the normalized similarity matrix - W = D^(-1/2) * A * D^(-1/2) */
double** create_W_mat(double **points, int n, int d){
    double **A_mat, **D_mat;
    double **first_multipication, **second_multipication;
    
    first_multipication = create_zero_mat(n, n);
    if (first_multipication == NULL)
        return NULL;

    second_multipication = create_zero_mat(n, n);
    if (second_multipication == NULL){
        freeMat(first_multipication, n);
        return NULL;
    }

    A_mat = create_A_mat(points, n, d);
    if (A_mat == NULL){
        freeMat(first_multipication, n);
        freeMat(second_multipication, n);
        return NULL;
    }

    D_mat = D_mat_normalized(points, n, d);
    if (D_mat == NULL){
        freeMat(first_multipication, n);
        freeMat(second_multipication, n);
        freeMat(A_mat, n);
        return NULL;
    }
    
    mat_multipication(D_mat, A_mat, first_multipication, n, n, n);
    mat_multipication(first_multipication, D_mat, second_multipication, n, n, n);

    /*free memory*/
    freeMat(A_mat, n);
    freeMat(D_mat, n);
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
    int i, j;
    
    for(i = 0; i < n; i++)
        for( j = 0; j < d; j++)
            matT[j][i] = mat[i][j];
}

/* return frobenius norm of matrice */
double frobenius_norm(double **matA,double **matB, int n, int k){
    double dist = 0.0;
    int i, j;
    
    for(i = 0; i < n; i++)
        for(j = 0; j < k; j++)
            dist += (matA[i][j]-matB[i][j])*(matA[i][j]-matB[i][j]);
    
    return dist;
}

/* calculation to find the new H matrice */
double **calcH(double** matH, double** matW, int n, int k){
    int i, j;
    double B = 0.5;
    double **newH,** matWH,**matHH,**matHHH,**matHt;
    
    newH = create_zero_mat(n, k);
    if (newH == NULL)
        return NULL;
    
    matWH = create_zero_mat(n, k);
    if (matWH == NULL){
        freeMat(newH, n);
        return NULL;
    }
    
    matHH = create_zero_mat(n, n);
    if (matHH == NULL){
        freeMat(matWH, n);
        freeMat(newH, n);
        return NULL;
    }
    
    matHHH = create_zero_mat(n, k);
    if (matHHH == NULL){
        freeMat(matHH, n);
        freeMat(matWH, n);
        freeMat(newH, n);
        return NULL;
    }

    matHt = create_zero_mat(k, n);
    if (matHt == NULL){
        freeMat(matHHH, n);
        freeMat(matHH, n);
        freeMat(matWH, n);
        freeMat(newH, n);
        return NULL;
    }
    
    mat_multipication(matW, matH, matWH, n, n, k);
    transpose(matH, matHt, n, k);
    mat_multipication(matH, matHt, matHH, n, k, n);
    mat_multipication(matHH, matH, matHHH, n, n, k);
   
    for(i = 0; i < n; i++)
        for(j = 0; j < k; j++)
            newH[i][j] = matH[i][j]*(1 - B + B*(matWH[i][j]/matHHH[i][j]));

    freeMat(matWH,n);
    freeMat(matHH,k);
    freeMat(matHHH,n);
    freeMat(matHt,k);
    
    return newH;
}

/* return the updated H matrix */
double** update_H_mat(double **matH, double **matW, int n, int k, double eps, int iter){
    int count = 0;
    double **newH, **copyH;
    
    copyH = matH;
    newH = calcH(copyH, matW, n, k);
    while(frobenius_norm(copyH,newH,n,k) >= eps && count < iter){
        copyH = newH;
        newH = calcH(copyH, matW, n, k);
        count++;
    }
    
    freeMat(copyH, n);
    return newH;
}

/* returns the number of rows and columns where res[0] contains the number of rows and res[1] contains the number of columns */
int* getDimensions(char* file_name){
    FILE *fileData = fopen(file_name, "r");
    int rownum = 0, colnum = 0, *res = malloc(2 * sizeof(int)), dim = 0;
    char c = 0;

    if (res == NULL) 
        return NULL;

    /* get number of rows and columns */
    while (1){
        c = fgetc(fileData); 
        if(c != EOF){
            if(c == ',')
                colnum++;
            else if (c == '\n'){
                rownum++;
                dim = colnum;
                colnum = 0;
            }
        }
	else
	    break;
    }
    dim++;
        
    res[0] = rownum;
    res[1] = dim;
    fclose(fileData);

    return res;
}

/* read and process the file to get the data as a points matrix */
double** getPoints(char* file_name, int n, int d){
    FILE *fileData = fopen(file_name, "r");
    int i, j;
    double num = 0, **points;
    
    points = (double**) calloc(n, sizeof(double*));
    if (points == NULL)
        return NULL;

    /* parse the data into points matrix */
    for(i = 0; i < n; i++){
        points[i] = (double*) calloc(d, sizeof(double)); 
        if (points[i] == NULL )
            return NULL;

        for(j = 0; j < d; j++){
            if(fscanf(fileData,"%lf",&num)!= EOF)
                points[i][j] = num;
            fgetc(fileData);
        }
    }

    fclose(fileData);
    return points;
}


int main(int argc, char *argv[])
{
    int n, d, *dims;
    char *file_name;
    double **points = NULL, **res = NULL;

    if (3 != argc){
        puts("An Error Has Occurred");
        return 1;
    }
    
    file_name = argv[2];
    
    /* get the dimensions and data from the file */
    dims = getDimensions(file_name);
    if (dims == NULL) {
        puts("An Error Has Occurred");
        return 1;
    }
    
    n = dims[0];
    d = dims[1];
    free(dims);

    points = getPoints(file_name, n, d);
    if (points == NULL) {
        puts("An Error Has Occurred");
        return 1;
    }

    /* run the algorithm based on the users input */
    if (strcmp(argv[1], "sym") == 0){
        res = create_A_mat(points, n, d);
        if (res == NULL)
            return 1;
            
        printm(res, n, n);
        freeMat(res, n);
    }

    if (strcmp(argv[1], "ddg") == 0){
        res = create_D_mat(points, n, d);
        if (res == NULL)
            return 1;
            
        printm(res, n, n);
        freeMat(res, n);
    }

    if (strcmp(argv[1], "norm") == 0){
        res = create_W_mat(points, n, d);
        if (res == NULL)
            return 1;
        
        printm(res, n, n);
        freeMat(res, n);
    }
    
    freeMat(points, n);
    return 0;
}