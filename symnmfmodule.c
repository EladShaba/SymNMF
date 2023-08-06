#define PY_SSIZE_T_CLEAN
#include <Python.h>       /* MUST include <Python.h>, this implies inclusion of the following standard headers:                        <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */
#include <math.h>         /* include <Python.h> has to be before any standard headers are included */
#include <stdio.h>
#include "symnmf.h"


PyObject* C_to_Py_mat(double **c_mat, int n, int d)
{
    /* cast matrice in C to numpy*/
    PyObject *py_matrice, *row;
    int i,j;
    py_matrice = PyList_New(n);/*list size n, each block is d : matrice nxd*/
    if (NULL == py_matrice){
        return NULL;
    }
    for (i = 0; i<n; i++)
    {
        row = PyList_New(d);
        if (NULL == row){
            return NULL;
        }
        for (j = 0; j < d; j++){
            PyList_SET_ITEM(row, j, Py_BuildValue("d", c_mat[i][j]));/*d-Convert a C double to a Python float*/
        }
        PyList_SetItem(py_matrice, i, Py_BuildValue("O", row));/*O-Pass a Python object untouched*/
    }
    return py_matrice;
}

double* Py_to_C_array(PyObject *py_data, int n, int d)
{   /* cast python array to C 1D array */
    PyObject *point, *index;
    double *C_data;
    int i,j;
    int flag = 0;
    C_data = (double*) calloc(n*d,sizeof(double));
    for ( i = 0; i < n; i++) {
        point = PyList_GetItem(py_data, i); /* Return the object at position index in the list pointed to by list */
        if (!PyList_Check(point)){
            /*if not list object, go to next element*/
            continue;
        }
        /* save py_data in C_data */
        for ( j = 0; j<d; j++) {
            index = PyList_GetItem(point, j);
            C_data[i*d + j] = PyFloat_AsDouble(index);/* i*d+j since we transform to 1D c array*/
            if (PyErr_Occurred() && C_data[i*d + j]  == -1.0){ /*if faliure than flag and break*/
                flag = 1;
                break;
            }
        }
        if (1 == flag){
            break;
        }
    }
    return C_data;
}

static PyObject* sym(PyObject *self, PyObject *args)
{
    PyObject *Py_points, *res;
    int n,d;
    double *C_arr;
    double **points, **A_mat;
    if (!PyArg_ParseTuple(args, "Oii", &Py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(Py_points)){
        return NULL;
    }
    C_arr = Py_to_C_array(Py_points, n, d);
    if (NULL == C_arr){
        return NULL;
    }
    points = arr_cast_mat(C_arr,n,d);
    if (NULL == points){
        return NULL;
    }
    A_mat = create_A_mat(points, n,d);
    if (NULL == A_mat){
        return NULL;
    }
    res = C_to_Py_mat(A_mat, n, n); /*cast A to python result*/
    /*free all memory*/
    free(C_arr); /*check if need to iterate over all array*/
    freeMat(A_mat, n);
    freeMat(points,n); /*two matrices to free*/
    return res;
}

static PyObject* ddg(PyObject *self, PyObject *args)
{
    PyObject *Py_points, *res;
    int n,d;
    double *C_arr;
    double **points, **A_mat, **D_mat;
    if (!PyArg_ParseTuple(args, "Oii", &Py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(Py_points)){
        return NULL;
    }
    C_arr = Py_to_C_array(Py_points, n, d);
    if (NULL == C_arr){
        return NULL;
    }
    points = arr_cast_mat(C_arr,n,d);
    if (NULL == points){
        return NULL;
    }
    A_mat = create_A_mat(points, n,d);
    if (NULL == A_mat){
        return NULL;
    }
    D_mat =create_D_mat(A_mat,n);
    if (NULL == D_mat){
        return NULL;
    }
    res = C_to_Py_mat(D_mat, n, n); /*cast D to python resutly*/
    /*free all memory*/
    free(C_arr); /*check if need to iterate over all array*/
    freeMat(A_mat, n);
    freeMat(D_mat, n);
    freeMat(points,n); /*three matrices to free*/
    return res;
}

static PyObject* norm(PyObject *self, PyObject *args)
{
    PyObject *Py_points, *res;
    int n,d;
    double *C_arr;
    double **points, **W_mat;
    if (!PyArg_ParseTuple(args, "Oii", &Py_points, &n, &d)){
        return NULL;
    }
    if (!PyList_Check(Py_points)){
        return NULL;
    }
    C_arr = Py_to_C_array(Py_points, n, d);
    if (NULL == C_arr){
        return NULL;
    }
    points = arr_cast_mat(C_arr,n,d);
    if (NULL == points){
        return NULL;
    }
    W_mat = create_W_mat(points, n,d);
    if (NULL == W_mat){
        return NULL;
    }
    res = C_to_Py_mat(W_mat, n, n); /*cast D to python resutly*/
    /*free all memory*/
    free(C_arr); /*check if need to iterate over all array*/
    freeMat(W_mat, n); 
    freeMat(points,n); /*two matrices to free*/
    return res;
}

static PyObject* updateH(PyObject *self, PyObject *args)
{
    PyObject *Py_H, *Py_W, *res;
    int n, k, iter;
    double *C_arr_H, *C_arr_W;
    double **H_mat, **W_mat, eps;
    if (!PyArg_ParseTuple(args, "OOiidi", &Py_H, &Py_W, &n, &k, &eps, &iter)){
        return NULL;
    }
    if (!PyList_Check(Py_H) || !PyList_Check(Py_W)){
        return NULL;
    }
    C_arr_H = Py_to_C_array(Py_H, n, k);
    C_arr_W = Py_to_C_array(Py_W, n, n);
    if (NULL == C_arr_H || NULL == C_arr_W){
        return NULL;
    }
    H_mat = arr_cast_mat(C_arr_H, n, k);
    W_mat = arr_cast_mat(C_arr_W, n, n);
    if (NULL == H_mat || NULL == W_mat){
        return NULL;
    }
    H_mat = update_H_mat(H_mat, W_mat, n, k, eps, iter);
    if (NULL == H_mat){
        return NULL;
    }
    res = C_to_Py_mat(H_mat, n, k); /*cast D to python resutly*/
    /*free all memory*/
    free(C_arr_H);
    free(C_arr_W); /*check if need to iterate over all array*/
    freeMat(W_mat, n); 
    freeMat(H_mat, n); /*two matrices to free*/
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
    "symnmf", 
    NULL, 
    -1,  
    symnmfMethods 
};

PyMODINIT_FUNC
PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
