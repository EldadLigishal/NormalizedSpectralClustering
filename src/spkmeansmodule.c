#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <spkmeans.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#define LINESIZE 1000


// need to edit the whole file

static void fit(PyObject *self, PyObject *args) {
    PyObject *_matrix;
    PyObject *line;
    char* goal;
    double **matrix;
    int n, d, i, j;

    if (!PyArg_ParseTuple(args, "Oiis", &n, &d, &goal, &_matrix)) {
        return NULL;
    }
    if (!PyList_Check(_matrix)){
        return NULL;
    }
    if(k>n){
        printf("Invalid Input! \n");
        return NULL;
    }

    matrix = createMat(n, d);
    assert(inputMat!=NULL);


    for (i = 0; i < n; i++){
        line = PyList_GetItem(_inputMat, i);
        for(j = 0 ; j < d ; j++){
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            matrix[i][j] = obj;
        }
    }

    if (strcmp(goal, "wam") == 0){
        printWAM(matrix, d, n);
    }
    if (strcmp(goal, "ddg") == 0){
        printDDG(matrix, d, n);
    }
    if (strcmp(goal, "lnorm") == 0){
        printLNORM(matrix, d, n);
    }
    if (strcmp(goal, "jacobi") == 0){
        printJacobi(matrix, d, n);
    }
    
    freeMemory(matrix,n);
}


static PyObject* fit_spk(PyObject *self, PyObject *args){
    PyObject *_inputMat;
    PyObject *_clusters;
    PyObject *line;
    PyObject *result;
    double epsilon, obj, **inputMat, **clusters, *array;
    int i, j, len;

    /* This parses the Python arguments into a int (i) variable named k ,
     *  int (i) variable named max_itr, double (d) variable named epsilon,
     *  int (i) variable named n, double** (O) variable named input_matrix
     *  double** (O) variable named clusters .
     *  gets k, 100, EPSILON, n, d, matrix, centroids.tolist()
     */
    if (!PyArg_ParseTuple(args, "iiiiiO", &k, &max_itr, &epsilon, &n, &d, &_inputMat, &_clusters)) {
        return NULL;
    }
    if (!PyList_Check(_inputMat) || !PyList_Check(_clusters)){
        return NULL;
    }
    if(k>n){
        printf("Invalid Input! \n");
        return NULL;
    }

    if (k == 0) {
        k = getEigengapHeuristic(array, len);
    }

    inputMat = createMat(n, d);
    assert(inputMat!=NULL);

    for (i = 0; i < n; i++){
        line = PyList_GetItem(_inputMat, i);
        for(j = 0 ; j < d ; j++){
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            inputMat[i][j] = obj;
        }
    }

    /*
     *  initialize centroids µ1, µ2, ... , µK
     */
    clusters = createMat(k,d);
    assert(clusters != NULL);

    for (i = 0; i < k; i++) {
        line = PyList_GetItem(_clusters, i);
        for(j=0 ; j<d ; j++) {
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            clusters[i][j] = obj;
        }
    }

    clusters = calculateCentroids(epsilon, inputMat, clusters);

    result = PyList_New(k);
    if(result == NULL){
        return NULL;
    }
    for(i = 0 ; i < k; i++){
        line = PyList_New(d);
        if(line==NULL){
            return NULL;
        }
        for(j = 0 ; j < d ; j++) {
            PyList_SetItem(line,j,PyFloat_FromDouble(clusters[i][j]));
        }
        PyList_SetItem(result, i, line);
    }

    freeMemory(inputMat,n);
    freeMemory(clusters,k);
    return Py_BuildValue("O", result);
}

static PyMethodDef myMethods[] = {
        { "fit",
        (PyCFunction)fit_spk, METH_VARARGS, PyDoc_STR("Input: Points, Centroids, Iterations and Clusters. Output: Centroids") },
        {"goal",
        (PyCFunction)fit, METH_VARARGS, PyDoc_STR("Input: Points, Centroids, Iterations and Clusters. Output: Centroids") },
        { NULL, NULL, 0, NULL }
};

static struct PyModuleDef myspkmeans = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans",
        NULL,
        -1,
        myMethods
};

PyMODINIT_FUNC PyInit_myspkmeans(void) {
    PyObject *m;
    m = PyModule_Create(&myspkmeans);
    if (!m) {
        return NULL;
    }
    return m;
}


/*
 *  creates 2-dimensional arrays
 */
double** createMat(int col, int row){
    int i;
    double ** matrix = (double**)malloc(col* sizeof(double *));
    assert(matrix != NULL);
    for(i=0; i < col; ++i){
        matrix[i]= (double*)malloc(row* sizeof(double));
        assert(matrix[i] != NULL);
    }
    return matrix;
}


double** calculateCentroids(double epsilon, double** inputMat, double** clusters){
    int i,j;
    /*
     *  groupOfClusters := group of clusters by S1,...,SK, each cluster Sj is represented by it’s
     *    centroid  which is the mean µj ∈ Rd of the cluster’s members.
     */
    double** groupOfClusters = NULL;
    /*
     *  groupOfClusters := [[0.0,0.0,0.0,0,0,0.0]
     *                     ,[0.0,0.0,0.0,0,0,0.0]]
     */
    groupOfClusters = createMat(k, n);
    assert(groupOfClusters != NULL);
    for(i=0; i<k; i++){
        for(j=0;j<n;j++){
            groupOfClusters[i][j] = 0.0;
        }
    }
    algorithm(clusters,inputMat,groupOfClusters, epsilon);

    /*
     * freeing memory
     */
    freeMemory(groupOfClusters, k);
    return clusters;
}

