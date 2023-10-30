#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>        
#include <stdio.h>
#include "symnmf.h"

/*cast a C matrice to python matrice*/
PyObject* C_to_Py_mat(double **c_mat, int n, int d){
    PyObject *py_matrice, *row;
    int i, j;
    
    if (c_mat == NULL)
        return NULL;

    /*parse a n x d C matrice to a n x d python matrice*/
    py_matrice = PyList_New(n); 
    if (py_matrice == NULL)
        return NULL;

    for (i = 0; i < n; i++){
        row = PyList_New(d);
        if (row == NULL)
            return NULL;

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
    if (C_mat == NULL)
        return NULL;

    for (i = 0; i < n; i++) {
        C_mat[i] = malloc(d * sizeof(double));
        if (C_mat[i] == NULL)
            return NULL;

        py_row = PyList_GetItem(py_data, i);
        for (j = 0; j < d; j++)
            C_mat[i][j] = PyFloat_AsDouble(PyList_GetItem(py_row, j));
    }

    return C_mat;
}

/* create and return a Python A matrix */
static PyObject* sym(PyObject *self, PyObject *args){
    PyObject *Py_points, *res;
    int n, d;
    double **points, **A_mat;
    
    if (!PyArg_ParseTuple(args, "O", &Py_points))
        return NULL;
    if (!PyList_Check(Py_points))
        return NULL;

    n = PyList_Size(Py_points);
    d = PyList_Size(PyList_GetItem(Py_points, 0));
    
    points = Py_to_C_mat(Py_points, n, d);
    if (points == NULL)
        return NULL;
    
    A_mat = create_A_mat(points, n, d);
    if (A_mat == NULL){
        freeMat(points, n);
        return NULL;
    }
    res = C_to_Py_mat(A_mat, n, n);
    
    /*free all memory*/
    freeMat(A_mat, n);
    freeMat(points, n);
    
    return res;
}

/* create and return a Python D matrix */
static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *Py_points, *res;
    int n, d;
    double **points, **D_mat;
    
    if (!PyArg_ParseTuple(args, "O", &Py_points))
        return NULL;
    if (!PyList_Check(Py_points))
        return NULL;
    
    n = PyList_Size(Py_points);
    d = PyList_Size(PyList_GetItem(Py_points, 0));

    points = Py_to_C_mat(Py_points, n, d);
    if (points == NULL)
        return NULL;
    
    D_mat = create_D_mat(points, n, d);
    if (D_mat == NULL){
        freeMat(points, n);
        return NULL;
    }

    res = C_to_Py_mat(D_mat, n, n);
    
    /*free all memory*/
    freeMat(D_mat, n);
    freeMat(points, n);
    
    return res;
}

/* create and return a Python W matrix */
static PyObject* norm(PyObject *self, PyObject *args){
    PyObject *Py_points, *res;
    int n, d;
    double **points, **W_mat;
    
    if (!PyArg_ParseTuple(args, "O", &Py_points))
        return NULL;
    if (!PyList_Check(Py_points))
        return NULL;
    
    n = PyList_Size(Py_points);
    d = PyList_Size(PyList_GetItem(Py_points, 0));

    points = Py_to_C_mat(Py_points, n, d);
    if (points == NULL)
        return NULL;
    
    W_mat = create_W_mat(points, n, d);
    if (W_mat == NULL){
        freeMat(points, n);
        return NULL;
    }

    res = C_to_Py_mat(W_mat, n, n);    
    
    /*free all memory*/
    freeMat(W_mat, n); 
    freeMat(points, n);
    
    return res;
}

/* return the updated Python H matrix */
static PyObject* updateH(PyObject *self, PyObject *args){
    PyObject *Py_H, *Py_W, *res;
    int n, k, iter;
    double **H_mat, **W_mat, eps;
    
    if (!PyArg_ParseTuple(args, "O|O|i|d|i", &Py_H, &Py_W, &k, &eps, &iter))
        return NULL;
    if (!PyList_Check(Py_H) || !PyList_Check(Py_W))
        return NULL;

    n = PyList_Size(Py_W);

    H_mat = Py_to_C_mat(Py_H, n, k);
    if (H_mat == NULL)
        return NULL;
    
    W_mat = Py_to_C_mat(Py_W, n, n);
    if (W_mat == NULL){
        freeMat(H_mat, n);
        return NULL;
    }

    H_mat = update_H_mat(H_mat, W_mat, n, k, eps, iter);
    res = C_to_Py_mat(H_mat, n, k);

    /*free all memory*/
    freeMat(W_mat, n); 
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
