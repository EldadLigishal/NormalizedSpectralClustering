#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include "spkmeans.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

static PyObject* fitope(PyObject *self, PyObject *args);
static PyObject* fitspk(PyObject *self, PyObject *args);
static PyObject* fitope_help(PyObject* _matrix, int d, int n, char* goal);
static PyObject* fitspk_help(int k , int maxItr, int d, PyObject* _matrix, PyObject* _centroids);
static PyObject* fitspkgetT(PyObject *self, PyObject *args);
static PyObject* fitspkgetT_help(PyObject* _matrix, int k, int d, int n);

/**
 *
 * @param self
 * @param args
 */
static PyObject* fitope(PyObject *self, PyObject *args){
    PyObject *_matrix;
    int d,n;
    char* goal;
    if (!PyArg_ParseTuple(args, "Oiis", &_matrix, &d, &n, &goal)){
        return NULL;
    }
    if (!PyList_Check(_matrix)){
        return NULL;
    }
    return Py_BuildValue("O", fitope_help(_matrix, d, n, goal));

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
    int k, maxItr,d;

    if (!PyArg_ParseTuple(args, "iiiOO", &k, &maxItr, &d, &_matrix, &_centroids)){
        return NULL;
    }
    if (!PyList_Check(_matrix) || !PyList_Check(_centroids)){
        return NULL;
    }
    return Py_BuildValue("O", fitspk_help(k, maxItr, d, _matrix, _centroids));
}

/**
 *
 * @param _matrix
 * @param d
 * @param n
 * @param goal
 * @return
 */
static PyObject* fitope_help(PyObject* _matrix, int d, int n, char* goal) {
    PyObject *line;
    PyObject *result;
    double **matrix;
    double ** toRetern;
    double obj;
    int i, j;
    double* eigenvalues;
    double** V;

    toRetern = NULL;
    /*
     * initialize input matrix.
     * converting Object to double
     */
    matrix = createMat(n, d);
    if (matrix == NULL) {
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

    if (strcmp(goal, "wam") == 0) {
        toRetern = getWeightedMatrix(matrix, d, n);
        freeMemory(matrix, n);
    }
    if (strcmp(goal, "ddg") == 0) {
        toRetern = getDiagonalDegreeMatrix(matrix, d,n);
        freeMemory(matrix, n);
    }
    if (strcmp(goal, "lnorm") == 0){
        toRetern = getLaplacianMatrix(matrix, d, n);
        freeMemory(matrix, n);
    }
    if (strcmp(goal, "jacobi") == 0){
        eigenvalues = (double *) malloc(n* sizeof (double));
        if(!eigenvalues){
            printf("Invalid Input!");
            exit(0);
        }
        for (i = 0; i < n; i++) {
            eigenvalues[i] = 0.0;
        }

        V = jacobiAlgorithm(matrix,n,eigenvalues);
        toRetern = concatenation(V,eigenvalues,n);

        freeMemory(V,n);
        free(eigenvalues);
    }

    /*
     * converting toReturn matrix from double to object
     */
    if(strcmp(goal, "jacobi") == 0){
        result = PyList_New(n + 1);
        if(result == NULL){
            return NULL;
        }
        for(i=0;i<(n + 1);i++){
            line = PyList_New(n);
            if(line==NULL){
                return NULL;
            }
            for(j=0;j<n;j++){
                PyList_SetItem(line,j,PyFloat_FromDouble(toRetern[i][j]));
            }
            PyList_SetItem(result, i, line);
        }
        freeMemory(toRetern, n + 1);
    }
    else{
        result = PyList_New(n);
        if(result == NULL){
            return NULL;
        }
        for(i=0;i<n;i++){
            line = PyList_New(n);
            if(line==NULL){
                return NULL;
            }
            for(j=0;j<n;j++){
                PyList_SetItem(line,j,PyFloat_FromDouble(toRetern[i][j]));
            }
            PyList_SetItem(result, i, line);
        }
        freeMemory(toRetern, n);
    }
    return result;
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
static PyObject* fitspk_help(int k , int maxItr, int d, PyObject* _matrix, PyObject* _centroids){
    PyObject  *result;
    PyObject *line;
    int i, j, n;
    double **matrix, **centroids_temp;
    double obj;
    /*
     * initialize input matrix.
     * converting Object to double
     */
    n = (int) PyObject_Length(_matrix);
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
     *  initialize centroids ??1, ??2, ... , ??K
     *  converting Object to double
     */
    centroids_temp = createMat(k, d);
    if(centroids_temp == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < d; i++) {
        line = PyList_GetItem(_centroids, i);
        for(j = 0 ; j < k ; j++) {
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            centroids_temp[i][j] = obj;
        }
    }

    kmeans(k,maxItr,d,n,matrix, centroids_temp);

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
    k = K;
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