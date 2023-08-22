#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>        
#include <stdio.h>
#include "symnmf.h"

/*cast a C matrice to python matrice*/
PyObject* C_to_Py_mat(double **c_mat, int n, int d){
    PyObject *py_matrice, *row;
    int i, j;
    
    if (c_mat == NULL){
        printf("problem with c mat\n");
        return NULL;
    }

    /*parse a n x d C matrice to a n x d python matrice*/
    py_matrice = PyList_New(n); 
    if (py_matrice == NULL){
        printf("problem in creating python list 1\n");
        return NULL;
    }

    for (i = 0; i < n; i++){
        row = PyList_New(d);
        if (row == NULL){
            printf("problem in creating python list 2\n");
            return NULL;
        }
        for (j = 0; j < d; j++)
            PyList_SetItem(row, j, PyFloat_FromDouble(c_mat[i][j]));
        PyList_SetItem(py_matrice, i, row);
    }

    return py_matrice;
}

/*cast a python matrice to C matrice*/
double** Py_to_C_mat(PyObject *py_data, int n, int d){   
    PyObject *py_row;
    double **C_mat;
    int i, j;

    /*parse a n x d python matrice to a n x d C matrice*/ 
    C_mat = malloc(n * sizeof(double*));
    if (C_mat == NULL){
        printf("problem in creating c list 1\n");
        return NULL;
    }

    for (i = 0; i < n; i++) {
        C_mat[i] = malloc(d * sizeof(double));
        if (C_mat[i] == NULL){
            printf("problem in creating c list 2\n");
            return NULL;
        }
        py_row = PyList_GetItem(py_data, i);
        for (j = 0; j < d; j++)
            C_mat[i][j] = PyFloat_AsDouble(PyList_GetItem(py_row, j));
    }

    return C_mat;
}

static PyObject* sym(PyObject *self, PyObject *args){
    PyObject *Py_points, *res;
    int n, d;
    double **points, **A_mat;
    
    if (!PyArg_ParseTuple(args, "O|i|i", &Py_points, &n, &d))
        return NULL;
    if (!PyList_Check(Py_points))
        return NULL;
    
    points = Py_to_C_mat(Py_points, n, d);
    A_mat = create_A_mat(points, n, d);
    res = C_to_Py_mat(A_mat, n, n);
    
    /*free all memory*/
    if (A_mat != NULL)
        freeMat(A_mat, n);
    if (points != NULL)
        freeMat(points,n); /*two matrices to free*/
    
    return res;
}

static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *Py_points, *res;
    int n,d;
    double **points, **A_mat, **D_mat;
    
    if (!PyArg_ParseTuple(args, "O|i|i", &Py_points, &n, &d))
        return NULL;
    if (!PyList_Check(Py_points))
        return NULL;
    
    points = Py_to_C_mat(Py_points, n, d);
    A_mat = create_A_mat(points, n,d);
    D_mat =create_D_mat(A_mat,n);
    res = C_to_Py_mat(D_mat, n, n);
    
    /*free all memory*/
    if (A_mat != NULL)
        freeMat(A_mat, n);
    if (D_mat != NULL)
        freeMat(D_mat, n);
    if (points != NULL)
        freeMat(points,n);
    
    return res;
}

static PyObject* norm(PyObject *self, PyObject *args){
    PyObject *Py_points, *res;
    int n,d;
    double **points, **W_mat;
    
    if (!PyArg_ParseTuple(args, "O|i|i", &Py_points, &n, &d))
        return NULL;
    if (!PyList_Check(Py_points))
        return NULL;
    
    points = Py_to_C_mat(Py_points, n, d);
    W_mat = create_W_mat(points, n,d);
    res = C_to_Py_mat(W_mat, n, n);    
    
    /*free all memory*/
    if (W_mat != NULL)
        freeMat(W_mat, n); 
    if (points != points)
        freeMat(points, n);
    
    return res;
}

static PyObject* updateH(PyObject *self, PyObject *args){
    PyObject *Py_H, *Py_W, *res;
    int n, k, iter;
    double **H_mat, **W_mat, eps;
    
    if (!PyArg_ParseTuple(args, "O|O|i|i|d|i", &Py_H, &Py_W, &n, &k, &eps, &iter))
        return NULL;
    if (!PyList_Check(Py_H) || !PyList_Check(Py_W))
        return NULL;

    H_mat = Py_to_C_mat(Py_H, n, k);
    W_mat = Py_to_C_mat(Py_W, n, n);
    H_mat = update_H_mat(H_mat, W_mat, n, k, eps, iter);
    res = C_to_Py_mat(H_mat, n, k);

    /*free all memory*/
    if (W_mat != NULL)
        freeMat(W_mat, n); 
    if (H_mat != NULL)
        freeMat(H_mat, n);
    
    return res;
}

static PyMethodDef symnmfMethods[] = {
    {"sym",                   
      sym,
      METH_VARARGS,         
      PyDoc_STR("sym")},
    {"ddg",                   
      ddg,
      METH_VARARGS,         
      PyDoc_STR("ddg")},
    {"norm",                   
      norm,
      METH_VARARGS,         
      PyDoc_STR("norm")},
    {"updateH",                   
      updateH,
      METH_VARARGS,         
      PyDoc_STR("updateH")}, 
    {NULL, NULL, 0, NULL}     
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "symnmfmod", 
    NULL, 
    -1,  
    symnmfMethods 
};

PyMODINIT_FUNC PyInit_symnmfmod(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
        return NULL;
    return m;
}
