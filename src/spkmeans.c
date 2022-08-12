#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <spkmeans.h>

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


/*
 *  create 2-dimensional arrays
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
                mat[i][j] = matrix1[i][k] * matrix2[k][j]
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