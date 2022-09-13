#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include "spkmeans.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

static void fitope(PyObject *self, PyObject *args);
static PyObject* fitspk(PyObject *self, PyObject *args);
static void fitope_help(PyObject* _matrix, int d, int n, char* goal);
static PyObject* fitspk_help(int k , int maxItr, int n, int d, PyObject* _matrix, PyObject* _centroids);
static PyObject* fitspkgetT(PyObject *self, PyObject *args);
static PyObject* fitspkgetT_help(PyObject* _matrix, int k, int d, int n);

/**
 *
 * @param self
 * @param args
 */
static void fitope(PyObject *self, PyObject *args){
    PyObject *_matrix;
    int d,n;
    char* goal;
    if (!PyArg_ParseTuple(args, "Oiis", &_matrix, &d, &n, &goal)){
        return;
    }
    if (!PyList_Check(_matrix)){
        return;
    }
    fitope_help(_matrix, d, n, goal);
}

/**
 *
 * @param self
 * @param args
 * @return
 * goal = spk
 */
static PyObject* fitspk(PyObject *self, PyObject *args) {
    PyObject *_matrix;
    PyObject *_centroids;
    int k, maxItr,d,n;
    
    if (!PyArg_ParseTuple(args, "iiiiOO", &k, &maxItr, &n, &d, &_matrix, &_centroids)){
        return NULL;
    }
    if (!PyList_Check(_matrix) || !PyList_Check(_centroids)){
        return NULL;
    }
    return Py_BuildValue("O", fitspk_help(k, maxItr, n, d, _matrix, _centroids));
}

/**
 *
 * @param _matrix
 * @param d
 * @param n
 * @param goal
 * @return
 */
static void fitope_help(PyObject* _matrix, int d, int n, char* goal){
    PyObject *line;
    double** matrix;
    double obj;
    int i, j;        

    /*
     * initialize input matrix.
     * converting Object to double
     */
    matrix = createMat(n, d);
    if (matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(0);
    }
    for (i = 0; i < n; i++) {
        line = PyList_GetItem(_matrix, i);
        for (j = 0; j < d; j++) {
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            matrix[i][j] = obj;
        }
    }

    if (strcmp(goal, "wam") == 0){
        printWAM(matrix,d,n);
    }
    if (strcmp(goal, "ddg") == 0){
        printDDG(matrix,d,n);
    }
    if (strcmp(goal, "lnorm") == 0){
        printLNORM(matrix,d,n);
    }
    if (strcmp(goal, "jacobi") == 0){
        printJACOBI(matrix,d,n);
    }
    /*
     * converting toReturn matrix from double to object
     */
    if(strcmp(goal, "jacobi") != 0){
        freeMemory(matrix, n);
    }
}


/**
 *
 * @param k
 * @param maxItr
 * @param n
 * @param d
 * @param _matrix
 * @param _centroids
 * @return
 */
static PyObject* fitspk_help(int k , int maxItr, int n, int d, PyObject* _matrix, PyObject* _centroids){
    PyObject  *result;
    PyObject *line;
    int i,j;
    double **matrix, **centroids;
    double obj;
    /*
     * initialize input matrix.
     * converting Object to double
     */
    matrix = createMat(n, d);
    if(matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        line = PyList_GetItem(_matrix, i);
        for (j = 0; j < d; j++) {
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            matrix[i][j] = obj;
        }
    }
    /*
     *  initialize centroids µ1, µ2, ... , µK
     *  converting Object to double
     */
    centroids = createMat(k, d);
    if(centroids == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < k; i++) {
        line = PyList_GetItem(_centroids, i);
        for(j = 0 ; j < d ; j++) {
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            centroids[i][j] = obj;
        }
    }

    controlPanel(k,maxItr,d,n,matrix, centroids);

    result = PyList_New(k);
    if(result == NULL){
        return NULL;
    }
    for(i=0;i<k;i++){
        line = PyList_New(d);
        if(line==NULL){
            return NULL;
        }
        for(j=0;j<d;j++){
            PyList_SetItem(line,j,PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_SetItem(result, i, line);
    }

    freeMemory(matrix,n);
    freeMemory(centroids,k);

    return result;
}


/**
 *
 * @param self
 * @param args
 * @return
 */
static PyObject* fitspkgetT(PyObject *self, PyObject *args) {
    PyObject *_matrix;
    int k, n, d;
    if (!PyArg_ParseTuple(args, "Oiii", &_matrix, &k, &d, &n)){
        return NULL;
    }
    if (!PyList_Check(_matrix)){
        return NULL;
    }
    return Py_BuildValue("O", fitspkgetT_help(_matrix,k,d,n));
}
/**
 *
 * @param _matrix
 * @param k
 * @param d
 * @param n
 * @return
 */
static PyObject* fitspkgetT_help(PyObject* _matrix, int k, int d, int n){
    PyObject *result;
    PyObject *line;
    double** matrix;
    double** tMatrix;
    double obj;
    int i,j;

    /*
     * initialize input matrix.
     * converting Object to double
     */
    matrix = createMat(n, d);
    if (matrix == NULL){
        printf("Invalid Input!");
        return 0;
    }
    for(i=0;i<n;i++){
        line = PyList_GetItem(_matrix, i);
        for (j=0;j<d;j++){
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            matrix[i][j] = obj;
        }
    }
    /*
     * calling getT function that are in spkmeans.c
     */
    tMatrix = getTMatrix(matrix,d,n,k);

    /*
     * TODO
     * we have to update the K !!!!!!!!!!!!!!!!!!!!!!!!
     */
    result = PyList_New(n);
    for(i=0;i<n;i++){
        line = PyList_New(k);
        for(j=0;j<k;j++)
        {
            PyList_SetItem(line,j,PyFloat_FromDouble(tMatrix[i][j]));
        }
        PyList_SetItem(result,i,line);
    }

    freeMemory(matrix,n);
    freeMemory(tMatrix,n);
    return result;
}


static PyMethodDef capiMethods[] = {
        { "fitope",
          (PyCFunction)fitope,METH_VARARGS,PyDoc_STR("centroids for k clusters")},
        { "fitspkgetT",
          (PyCFunction)fitspkgetT,METH_VARARGS,PyDoc_STR("centroids for k clusters")},
        { "fitspk",
          (PyCFunction)fitspk,METH_VARARGS,PyDoc_STR("centroids for k clusters")},
        {NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans",
        NULL,
        -1,
        capiMethods
};

PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject *m;
    m=PyModule_Create(&moduledef);
    if(!m){
        return NULL;
    }
    return m;
}