#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <spkmeans.h>

#define EPSILON pow(10, -5)

enum GOAL {
    SPK,
    WAM,
    DDG,
    LNORM,
    JACOBI
};


/*
 * argc := number of inputs.
 * argv := ["...",k,goal,input_filename].
 * argv[0] is the name of the program.
 */
int main(int argc, char* argv[]) {
    FILE *input_file;
    int k, i, n, j;
    char* goal;
    enum GOAL g;

    if (argc != 3) {
        printf("Invalid Input!");
        return 0;
    }

    k = atoi(argv[1]);
    goal = argv[2];
    input_file = fopen(argv[3], "r");

    if (input_file == NULL) {
        printf("Invalid Input!");
        return 0;
    }
    
    // good approach?
    if (strcmp(goal, "spk") == 0) {
        g = SPK;
    }

    if (strcmp(goal, "wam") == 0) {
        g = WAM;
    }

    if (strcmp(goal, "ddg") == 0) {
        g = DDG;
    }

    if (strcmp(goal, "lnorm") == 0) {
        g = LNORM;
    }

    if (strcmp(goal, "jacobi") == 0) {
        g = JACOBI;
    }

    if ((g != SPK) && (g != WAM) && (g != DDG) && 
        (g != LNORM) && (g != JACOBI)) {
        printf("Invalid Input!");
        return 0;
    }
    return 0;
}


double weight(double** matrix, int index1, int index2, int dim) {
    int i;
    double value = 0;
    double distance;
    for (i = 0; i < dim; i++) {
        distance = fabs(matrix[index1][i] - matrix[index2][i]);
        value = pow(distance, 2.0);
    }
    value = exp(-((sqrt(value)) / 2));
    return value;
}


/*
 * creates weighted matrix
 */
double** getWeightedMatrix(double** matrix, int cols, int rows) {
    int i, j;
    double** wMatrix;
    double value;
    wMatrix = createMat(cols, rows);
    for (i = 0; i < rows; i++) {
        for (j = 0; j < rows; j++) {
            value = 0;
            if (i != j) {
                value = weight(matrix,i , j, cols);
            }
            wMatrix[i][j] = value;
        }
    }
    return wMatrix;
}


/*
 *  creates diagonal matrix
 */
double** getDiagonalMatrix(double** matrix, int dim) {
    int i, j;
    double** dMatrix = NULL;
    dMatrix = createMat(dim, dim);
    for (i = 0; i < dim; ++i) {
        double sum = 0;
        for (j = 0; j < dim; j++) {
            sum += matrix[i][j];
        }
        dMatrix[i][i] = sum;
    }
    return dMatrix;
}

/*
 *  creates diagonal degree matrix
 */
double** getDiagonalD(double** wMatrix, int dim) {
    int i, j;
    double** dmatrix = NULL;
    dmatrix = getDiagonalMatrix(wMatrix, dim);
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; j++) {
            if (dmatrix[i][j] != 0) {
                dmatrix[i][j] = pow(dmatrix[i][j], -0.5);
            }
        }
    }
    return dmatrix;
}


/*
 *  creates normalized laplacian matrix
 */
double** getLaplacianMat(double** weightMatrix, double** diagMatrix, int dim) {
    int i, j;
    double** temp = multiplyMatrices(diagMatrix, weight, dim);
    double** lMatrix = multiplyMatrices(temp, diagMatrix, dim);
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (i == j) {
                lMatrix[i][j] = 1 - lMatrix[i][j];
            } else {
                lMatrix[i][j] = -(lMatrix[i][j]);
            }
        }
    }
    freeMemory(temp, dim);
    return lMatrix;
}


/*
 *  Part of step 5
 *  calculate off(A)^2
 */
double offOfMat(double** matrix, int dim) {
    double value = 0;
    int i, j;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (i != j) {
               value += pow(matrix[i][j], 2.0);
            }
        }
    }
    return value;
}


/*
 *  Jacobi algorithm
 */
double** jacobiAlgorithm(double** matrix, int dim) {
    double offset = 0;
    int iter = 0;
    // convergence condition
    while (iter < 0 || EPSILON < offset ) {

    }
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


double** transpose(double** inputMat, int dim) {
    int i, j;
    double** mat;
    mat = createMat(dim, dim);
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            mat[i][j] = inputMat[j][i];
        }
    }
    return mat;
}


/*
 *  multiply 2 matrices
 */
double** multiplyMatrices(double** matrix1, double** matrix2, int dim) {
    int j, i, k;
    double** mat;
    mat = createMat(dim, dim);
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            for (k = 0; k < dim; ++k) {
                mat[i][j] = matrix1[i][k] * matrix2[k][j];
          }
        }
    }
    return mat;
}


/*
 *  free memory function
 */
void freeMemory(double** matrix ,int len){
    int i;
    if(matrix == NULL){
        return;
    }
    for(i = 0; i < len ; i++){
        if(matrix[i] == NULL){
            continue;
        }
        free(matrix[i]);
    }
    free(matrix);
}