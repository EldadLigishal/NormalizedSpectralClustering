#include <Python.h>
#include <stdio.h>
#define LINESIZE 1000
#define PY_SSIZE_T_CLEAN


// need to edit the whole file


static PyObject* fit(PyObject *self, PyObject *args){
    PyObject *_inputMat;
    PyObject *_clusters;
    PyObject *line;
    PyObject *result;
    double obj;
    double **inputMat;
    double **clusters;
    double epsilon;
    int i;
    int j;

    /* This parses the Python arguments into a int (i) variable named k ,
     *  int (i) variable named max_itr, double (d) variable named epsilon,
     *  int (i) variable named n, double** (O) variable named input_matrix
     *  double** (O) variable named clusters .
     */
    if (!PyArg_ParseTuple(args, "iidiiOO", &k, &max_itr, &epsilon, &n, &d, &_inputMat, &_clusters)) {
        return NULL;
    }
    if (!PyList_Check(_inputMat) || !PyList_Check(_clusters)){
        return NULL;
    }
    if(k>n){
        printf("Invalid Input! \n");
        return NULL;
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
        (PyCFunction)fit, METH_VARARGS, PyDoc_STR("Input: Points, Centroids, Iterations and Clusters. Output: Centroids") },
        { NULL, NULL, 0, NULL }
};

static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        "kmeans module",
        -1,
        myMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
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