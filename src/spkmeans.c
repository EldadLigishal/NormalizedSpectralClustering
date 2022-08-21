#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "assert.h"
#include <string.h>
#include "spkmeans.h"


/*
 * argc := number of inputs.
 * argv := ["...",k,goal,input_filename].
 * argv[0] is the name of the program.
 */
int main(int argc, char* argv[]) {
    FILE *input_file;
    int k, i, n, d, j;
    char* goal;
    GOAL g;
    double **inputMatrix;
    int validGoal =0;

    assert(argc > 0);
    if (argc != 4) {
        printf("Invalid Input!");
        return 0;
    }

    // read the input
    k = atoi(argv[1]);
    if (k < 0){
        printf("Invalid Input!");
        return 0;
    }
    goal = argv[2];
    input_file = fopen(argv[3], "r");
    if (input_file == NULL){
        printf("Invalid Input!");
        return 0;
    }

    /*
     * calculating the columns and rows of the input file
     * check if k==n?
     */
    n = calculateRows(input_file);  // columns
    d = calculateCol(input_file);   // rows

    /*
     * build a matrix from input file
     */
    inputMatrix = createMat(n, d);
    assert(inputMatrix != NULL);
    fclose(input_file);

    double **weightedMatrix = getWeightedMatrix(inputMatrix, n, d);
    double **diagMatrix = getDiagonalD(weightedMatrix, d);
    double **lapMatrix = getLaplacianMat(weightedMatrix, diagMatrix, d);

    if (strcmp(goal, "spk") == 0){
        g = SPK;
        validGoal=1;
    }
    if (strcmp(goal, "wam") == 0){
        g = WAM;
        validGoal=1;
    }
    if (strcmp(goal, "ddg") == 0){
        g = DDG;
        validGoal=1;
    }
    if (strcmp(goal, "lnorm") == 0){
        g = LNORM;
        validGoal=1;
    }
    if (strcmp(goal, "jacobi") == 0){
        g = JACOBI;
        validGoal=1;
    }
    if(validGoal == 0){
        printf("Invalid Input!");
        return 0;
    }

    freeMemory(inputMatrix, n);
    return 0;
}

/*
 * calculate the value of exp(−||xi − xj||/2)
 */
double calculateWeight(double** matrix, int index1, int index2, int dim){
    int i;
    double value = 0;
    double distance;
    for (i = 0; i < dim; i++) {
        distance = fabs(matrix[index1][i] - matrix[index2][i]);
        value += pow(distance, 2.0);
    }
    value = exp((-1)*( (sqrt(value)) / 2));
    return value;
}

/*
 * Form the weighted adjacency matrix W ∈ R^(n×n).
 * The weights are symmetric (wij = wji) and non-negative (wij ≥ 0).
 * wii = 0 for all i’s
 * the rest of the values are set to: wij = exp(−||xi − xj||/2)
 */
double** getWeightedMatrix(double** matrix, int n){
    int i, j;
    double** wMatrix;
    wMatrix = createMat(n, n);
    if(!wMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    for (i = 0; i < n; i++) {
        for (j = (i+1); j < n; j++) {
            wMatrix[i][j] = calculateWeight(matrix,i , j, n);
            wMatrix[j][i] = wMatrix[i][j];
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
 * check if a matrix is diagonal
 */
int isDiagonalMatrix(double** matrix, int dim) {
    int i, j;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if ((i != j) && matrix[i][j] != 0) {
               return 0;    //return false
            }
        }
    }
    return 1;   //return true
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

int getSign(double number) {
    if (number >= 0) {
        return 1;
    }
    return -1;
}


double** getUnitMatrix(int dim) {
    int i, j;
    double **matrix = createMat(dim, dim);
    for ( i = 0; i < dim; i++) {
        for ( j = 0; j < dim; i++) {
            if (i == j) {
                matrix[i][j] = 1;
            } else {        // TODO: CHANGE CREATEMATRIX TO CALLOC
                matrix[i][j] = 0;
            }
        }
    }
    return matrix;
}


double** getRotationMat(double **matrix, int dim) {
    int i, j;
    int index_i = 0;
    int index_j = 1;
    int sign_theta;
    double t, c, s;
    double **pMat;
    double theta = 0;
    double maxValue = matrix[0][1];
    
    for (i = 0; i < dim - 1; i++) {
        for ( j = i + 1; i < dim; j++) {
            if (maxValue < matrix[i][j]) {
                maxValue = matrix[i][j];
                index_i = i;
                index_j = j;
            }
        }
    }

    theta = ((matrix[index_j][index_j] - matrix[index_i][index_i])
     / 2*matrix[index_i][index_j]);
    sign_theta = getSign(theta);
    
    t = (sign_theta / (fabs(theta) + sqrt(pow(theta, 2) + 1)));
    c = (1 / sqrt(pow(t, 2) + 1));
    s = c * t;
    pMat = getUnitMatrix(dim);
    pMat[index_i][index_i] = c;
    pMat[index_i][index_j] = s;
    pMat[index_j][index_j] = c;
    pMat[index_j][index_i] = -s;

    return pMat;
}

/*
 *  Jacobi algorithm
 */
double** jacobiAlgorithm(double** matrix, int dim) {
    double offset = 0;
    int iter = 0;


    // building a rotation matrix P
    double **pMatrix = getRotaionMat(matrix, dim);      //TODO: Check if we are in even iteration?
    double **vMatrix = getRotaionMat(matrix, dim);
    // building matrix A'
    double **aTagTemp = multiplyMatrices(transpose(pMatrix, dim), matrix, dim);
    double **aTagMat = multiplyMatrices(aTagTemp, pMatrix, dim);

    offset = offOfMat(matrix, dim) - offOfMat(aTagMat, dim);

    // repeat until A' is diagonal matrix
    while ((isDiagonalMatrix(aTagMat, dim) == 0) || (iter < 100 || EPSILON < offset )) {
        pMatrix = getRotaionMat(aTagMat, dim);
        vMatrix = multiplyMatrices(vMatrix, pMatrix, dim);
        aTagTemp = multiplyMatrices(transpose(pMatrix, dim), aTagMat, dim);
        aTagMat = multiplyMatrices(aTagTemp, pMatrix, dim);
        offset = offOfMat(matrix, dim) - offOfMat(aTagMat, dim);
    }

    double **aMatrix = aTagMat;


    freeMemory(aTagTemp, dim);
}


/*
 *  creates 2-dimensional arrays
 *  NEEDS TO CHANGE TO CALLOC
 *  TODO matrix full of zeroz
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
                mat[i][j] += matrix1[i][k] * matrix2[k][j];
          }
        }
    }
    return mat;
}


/*
 * calculation the number of col in fileName
 */
int calculateCol(char* fileName){
    int cnt = 0;
    FILE* ifp;
    char* token;
    const char breaks[] = ",";
    char line[LINESIZE];
    /*
     * open file
     */
    ifp = fopen(fileName,"r");
    if(ifp == NULL) {
        printf("Invalid Input! \n");
        return -1;
    }
    /*
     * calculating the number of col
     */
    fgets(line, 1000, ifp);
    token = strtok(line, breaks);
    while(token != NULL)
    {
        token = strtok(NULL, breaks);
        ++cnt;
    }

    /*
     * close file
     */
    fclose(ifp);
    return cnt;

}


/*
 * calculate the number of rows in fileName
 */
int calculateRows(char* fileName){
    int cnt = 0;
    char line[LINESIZE];
    FILE* ifp;

    /*
     * open file
     */
    ifp = fopen(fileName,"r");
    if(ifp ==NULL) {
        printf("Invalid Input! \n");
        return 0;
    }
    /*
     * calculating the number of rows
     */
    while (fgets(line, LINESIZE, ifp) != NULL){
        cnt++;
    }
    /*
     * close file
     */
    fclose(ifp);
    return cnt ;
}


/*
 * fill 2-dimensional arrays
 */
void fillMat(char* fileName,double** inputMat){
    FILE* ifp;
    char* token;
    const char breaks[] = ",";
    char line[LINESIZE];
    int row=0,col=0;
    char* useless=NULL;

    ifp = fopen(fileName,"r");
    if(ifp == NULL) {
        printf("Invalid Input! \n");
        return;
    }
    while (fgets(line,LINESIZE,ifp) != NULL){
        token = strtok(line,breaks);
        while (token != NULL){
            inputMat[row][col] = strtod(token, &useless);
            token = strtok(NULL,breaks);
            col++;
        }
        col=0;
        row++;
    }
    fclose(ifp);
}


void printWAM(double** matrix, int col,int row) {
    double **wMatrix = getWeightedMatrix(matrix, col, row);
    printMat(wMatrix, row, row);
    freeMemory(wMatrix, row);
}


void printDDG(double** matrix, int col,int row) {
    double **wMatrix = getWeightedMatrix(matrix, col, row);
    double **dMatrix = getDiagonalD(matrix, row);
    printMat(dMatrix, row, row);
    freeMemory(wMatrix, row);
    freeMemory(dMatrix, row);
}


void printLNORM(double** matrix, int col, int row) {
    double **wMatrix = getWeightedMatrix(matrix, col, row);
    double **dMatrix = getDiagonalD(matrix, row);
    double **lMatrix = getLaplacianMat(wMatrix, dMatrix, row);
    printMat(lMatrix, row, row);
    freeMemory(wMatrix, row);
    freeMemory(dMatrix, row);
    freeMemory(lMatrix, row);
}

void printJACOBI(double** matrix, int col,int row) {}


void printMat(double** matrix, int col,int row){
    int i, j;
    
    printf("\nPrinting the 2D Array\n");
    for(i = 0; i < j; i++)
    {
        for(j = 0; j < i; j++)
        {
            printf("%.4f|", matrix[i][j]);
        }
        printf("\n\n");
    }
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
