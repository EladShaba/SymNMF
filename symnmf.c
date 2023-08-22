#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "symnmf.h"

/*create a n x d zero matrice*/
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

/*free memory of matrice*/
void freeMat(double **mat, int n){ 
    int i;
    for(i = 0; i < n; i++)
        free(mat[i]);
    free(mat);
}

/*create the A-Simlarity matrice*/
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

/*create the D-Diagonal Degree Matrix*/
double** create_D_mat(double **A_mat, int n){
    double **D_mat, sum;
    int i ,j;

    if (A_mat == NULL)
        return NULL;

    D_mat = create_zero_mat(n, n);
    if (D_mat == NULL)
        return NULL;

    for (i = 0; i < n; i++){
        sum = 0.0;
        for (j = 0; j < n; j++)
            sum += A_mat[i][j];
        D_mat[i][i] = sum;
    }

    return D_mat;
}

/*create the normalized D-Diagonal Degree Matrix (D^(-1/2))*/
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

/*fill C with the result of A*B where
A is a n x m matrice, B is a m x d matrice and C is a n x d matrice*/ 
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

double** create_W_mat(double **points, int n, int d){
    /* return W - the normalized similarity matrix - W = D^(-1/2) * A * D^(-1/2) */
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

/*prints a n x d matrice*/
void printm(double **mat, int n, int d){
    int i,j;
    
    for ( i = 0; i < n; i++){
        for ( j = 0; j < d - 1; j++)
            printf("%.4f,", mat[i][j]);
        printf("%.4f\n", mat[i][d-1]);
    }
}

/*make matT be the transposed matrice of mat (which is a n x d matrice)*/
void transpose(double **mat, double **matT, int n, int d){
    int i,j;
    
    for(i = 0; i < n; i++)
        for( j = 0; j < d; j++)
            matT[j][i] = mat[i][j];
}

/*return frobenius norm of matrice*/
double frobenius_norm(double **mat, int n, int k){
    int i, j;
    double norm = 0;
    
    for(i = 0; i < n; i++)
        for(j = 0; j < k; j++)
            norm += pow(mat[i][j], 2);
    
    return norm;
}

/*creates a n x d matrice from the points in the file*/
double** get_points_input(FILE *ifp, int n, int d){
    int i, j;
    double **mat, point;

    mat = create_zero_mat(n,d);
    if (mat == NULL)
        return NULL;

    for(i = 0; i < n; i++)
        for(j = 0; j < d; j++)
            if (fscanf(ifp, "%lf", &point) != EOF)
                mat[i][j] = point;
    
    return mat;
}

double** update_H_mat(double **matH, double **matW, int n, int k, double eps, int iter){
    double **new_H, **numerator_H, **denominator_H, **first_mult, **mat_H_T, **frobenius_mat, norm;
    int i, j;
    
    if (matH == NULL || matW == NULL)
        return NULL;
    
    new_H = create_zero_mat(n, k);
    frobenius_mat = create_zero_mat(n, k);
    numerator_H = create_zero_mat(n, k);
    first_mult = create_zero_mat(k, k);
    denominator_H = create_zero_mat(n, k);
    mat_H_T = create_zero_mat(k, n);

    /*at least one of the matrices failed to load correctly*/
    if (numerator_H == NULL || mat_H_T == NULL || denominator_H == NULL
        || first_mult == NULL || new_H == NULL || frobenius_mat == NULL){
        if (numerator_H != NULL)
            freeMat(numerator_H, n);
        if (mat_H_T != NULL)
            freeMat(mat_H_T, n);
        if (first_mult != NULL)
            freeMat(first_mult, n);
        if (denominator_H != NULL)
            freeMat(denominator_H, n);
        if (new_H != NULL)
            freeMat(new_H, n);
        if (frobenius_mat != NULL)
            freeMat(frobenius_mat, n);
        return NULL;
    }


    while (1){
        transpose(matH, mat_H_T, n, k);
        mat_multipication(matW, matH, numerator_H, n, n, k);
        mat_multipication(mat_H_T, matH, first_mult, k, n, k);
        mat_multipication(matH, first_mult, denominator_H, n, k, k);

        /*parse the new H and parse the frobenius mat (newH - oldH)*/
        for (i = 0; i < n; i++)
            for (j = 0; j < k; j++){
                new_H[i][j] = matH[i][j] * (0.5 + ((0.5*numerator_H[i][j]) / denominator_H[i][j]));
                frobenius_mat[i][j] = new_H[i][j] - matH[i][j];
            }
        norm = frobenius_norm(frobenius_mat, n, k);
        
        if (iter == 0 || norm < eps)
            break;
        
        /*set matH to be the new_H*/
        freeMat(matH, n);
        matH = new_H;
        new_H = create_zero_mat(n, k);
        if (new_H == NULL)
            return NULL;
        
        iter -= 1;
    }
    
    /*free the matrices used for the operations*/
    freeMat(numerator_H, n);
    
    /*why isn't this working when uncommenting these lines????????*/
    //freeMat(mat_H_T, n);
    //freeMat(first_mult, n);
    
    freeMat(denominator_H, n);
    freeMat(frobenius_mat, n);
    freeMat(new_H, n);
    
    return matH;
}


int main(int argc, char** argv){
    /* if there are command-line arguments, they are interpered as filenames, and processed in order */
    FILE *ifp; 
    int n = 0, d = 1;   
    char *filename, charCount;
    double **points, **A_mat, **D_mat, **W_mat;

    filename = argv[2];
    ifp = fopen(filename, "r");
    if (ifp == NULL){
        printf("An Error Has Occurred");
        return 1;
    }

    while ((charCount = getchar()) != EOF){
        if (n == 0 && charCount == ',')
            d += 1;
        if (charCount == '\n')
            n += 1;
    }
    
    rewind(stdin); /* reset pointer*/

    /*load data-points from file into points n x d matrice*/
    points = get_points_input(ifp, n, d);
    if (points == NULL)
        return 1;

    fclose(ifp);

    if (strcmp(argv[1], "sym") == 0){
        A_mat = create_A_mat(points, n, d);
        if (A_mat == NULL)
            return 1;
        
        printm(A_mat, n, n);
        freeMat(A_mat, n);
    }

    if (strcmp(argv[1], "ddg") == 0){
        D_mat = create_D_mat(points, n);
        if (D_mat == NULL)
            return 1;
        
        printm(D_mat, n, n);
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
