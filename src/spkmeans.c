#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <string.h>
#include "spkmeans.h"
#include <stdbool.h>
#include "assert.h"

/*
 * argc := number of inputs.
 * argv := ["...", goal,filename].
 * argv[0] is the name of the program.
 * goal (enum) , file name (.txt or .csv)
 */
int main(int argc, char* argv[]) {
    char* filename;
    FILE* ifp;
    int n, d;
    char* goal;
    Goal g;
    double **inputMatrix;
    int validGoal =0;

    if (argc != 3) {
        printf("Invalid Input!");
        return 0;
    }
    /*
    * open file
    */
    filename = argv[2];
    ifp = fopen(filename, "r");
    if(ifp == NULL){
        printf("Invalid Input! \n");
        return 0;
    }

    goal = argv[1];
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
    /*
     * calculating the columns and rows of the input file
     * check if d==n?
     */
    n = calculateRows(filename);
    d = calculateCol(filename);
    /*
     * build a matrix from input file
     * starting by creating an empty matrix, and then we fill it.
     */
    inputMatrix = createMat(n, d);
    if (inputMatrix == NULL){
        printf("Invalid Input!");
        return 0;
    }
    fillMat(filename, inputMatrix);

    operation(inputMatrix,d,n,g);
    fclose(ifp);
    freeMemory(inputMatrix, n);
    return 0;
}


void operation(double **matrix, int dim,int num, Goal g) {
    if (g == WAM) {
        printWAM(matrix, dim,num);
    }
    if (g == DDG) {
        printDDG(matrix, dim, num);
    }
    if (g == LNORM) {
        printLNORM(matrix, dim, num);
    }
    if (g == JACOBI) {
        printJACOBI(matrix, dim);
    }
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
double** getWeightedMatrix(double** matrix, int dim, int num){
    int i, j;
    double** wMatrix;
    if(!matrix){
        printf("An Error Has Occured");
        exit(0);
    }
    wMatrix = createMat(num, num);
    if(!wMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    for (i = 0; i < num; i++) {
        for (j = (i+1); j < num; j++) {
            wMatrix[i][j] = calculateWeight(matrix,i , j, dim);
            wMatrix[j][i] = wMatrix[i][j];
        }
    }
    return wMatrix;
}

/*
 * Form the diagonal degree matrix D ∈ R^(n×n).
 * the diagonal equals to the sum of the i-th row of W.
 * if i=j dij = sumOf(wiz) for all z=1..n, otherwise dij=0
 */
double** getDiagonalDegreeMatrix(double** matrix,int dim, int num){
    int i, j;
    double sum;
    double** dMatrix;
    double** wMatrix;
    wMatrix = getWeightedMatrix(matrix, dim, num);
    dMatrix = createMat(num, num);
    if(!dMatrix){
        printf("An Error Has Occurred");
        exit(0);
    }
    for (i = 0; i < num; i++){
        sum = 0;
        for (j = 0; j < num; j++) {
            sum += wMatrix[i][j];
        }
        dMatrix[i][i] = sum;
    }
    freeMemory(wMatrix,num);
    return dMatrix;
}

/*
 * D^(1/2)
 * D is a Diagonal degree matrix
 */
double** getDiagonalMatrixPow2(double** matrix, int dim, int num){
    int i;
    double** dMatrix;
    dMatrix = getDiagonalDegreeMatrix(matrix, dim, num);
    for (i = 0; i < num; ++i) {
        dMatrix[i][i] = pow(sqrt(dMatrix[i][i]), -0.5);
    }
    return dMatrix;
}

/*
 * Form The Normalized Graph Laplacian matrix Lnorm ∈ R^(n×n).
 * Lnorm = I − (D^(-1/2) * W * D^(-1/2))
 *       = ( (Diagonal Degree Matrix * Weighted Matrix) * Diagonal Degree Matrix )
 */
double** getLaplacianMatrix(double** matrix, int dim, int num){
    int i, j;
    double** wMatrix;
    double** dMatrix;
    double** lnormMatrix;
    double** temp;
    wMatrix = getWeightedMatrix(matrix,dim, num);
    dMatrix = getDiagonalMatrixPow2(matrix,dim, num);
    /*
     * temp = D^(-1/2) * W
     */
    temp = multiplyMatrices(dMatrix, wMatrix, num);
    if(!temp){
        printf("An Error Has Occured");
        exit(0);
    }
    /*
     * lnormMatrix = (D^(-1/2) * W)*D^(-1/2)
     */
    lnormMatrix = multiplyMatrices(temp, dMatrix, num);
    if(!lnormMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    for (i = 0; i < num; i++) {
        for (j = 0; j < num; j++) {
            if (i == j) {
                lnormMatrix[i][j] = 1 - lnormMatrix[i][j];
            } else {
                lnormMatrix[i][j] = ((-1)*(lnormMatrix[i][j]));
            }
        }
    }
    freeMemory(temp, num);
    freeMemory(wMatrix,num);
    freeMemory(dMatrix,num);
    return lnormMatrix;
}
/*
 * check if a matrix is diagonal
 * diagonal matrix is a square matrix that consists of all zeros off the main diagonal.
 */
bool isDiagonalMatrix(double** matrix, int dim){
    int i, j;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if ((i != j) && matrix[i][j] != 0) {
               return false;
            }
        }
    }
    return true;
}

/*
 * sign(x) = 1 if x>=0, else 0
 */
int getSign(double number) {
    if (number >= 0) {
        return 1;
    }
    return -1;
}
/*
 * The Eigengap Heuristic
 * In order to determine the number of clusters k, we will use eigengap heuristic.
 */
int getEigengapHeuristic(double* array,int len){
    int maxIndex=0;
    int i;
    double deltaMax=0;
    /*
    * δ_max >= δi (= |λi − λi+1|)
    */
    for(i = 0; i < floor(len/2) ; i++){
        if( fabs(array[i]-array[i+1]) > deltaMax){
            maxIndex=i;
            deltaMax=fabs(array[i]-array[i+1]);
        }
    }
    return maxIndex+1;
}

double** getUnitMatrix(int dim) {
    int i, j;
    double **matrix = createMat(dim, dim);
    if(!matrix){
        printf("Invalid Input!");
        exit(0);
    }
    for ( i = 0; i < dim; i++) {
        for ( j = 0; j < dim; j++) {
            if (i == j) {
                matrix[i][j] = 1;
            } else {
                matrix[i][j] = 0;
            }
        }
    }
    return matrix;
}
/*
 *  creates 2-dimensional arrays
 *  NEEDS TO CHANGE TO CALLOC
 *  matrix full of zeros
 */
double** createMat(int col, int row){
    int i,j;
    double ** matrix = (double**)malloc(col* sizeof(double *));
    assert(matrix != NULL);
    for(i=0; i < col; ++i){
        matrix[i]= (double*)malloc(row* sizeof(double));
        assert(matrix[i] != NULL);
    }
    if(!matrix){
        printf("An Error Has Occurred");
        exit(0);
    }
    return matrix;
}
/*
 * Check if we have freed the input matrix
 */
double** transpose(double** matrix, int dim){
    int i, j;
    double** tMatrix;
    tMatrix = createMat(dim, dim);
    if(!tMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    for ( i = 0; i < dim; ++i)
        for ( j = 0; j < dim; ++j) {
            tMatrix[j][i] = matrix[i][j];
        }
    return tMatrix;
}

bool isSymmetric(double** matrix,int dim){
    /*
     * TODO
     * check if a matrix is symmetric or not
     */
}
/*
 *  multiplying matrices dimXdim
 *  the number of columns in the first matrix must be equal
 *  to the number of rows in the second matrix.
 */
double** multiplyMatrices(double** matrix1, double** matrix2, int dim) {
    int j, i, l;
    double** mat;
    mat = createMat(dim, dim);
    if(!mat){
        printf("An Error Has Occured");
        exit(0);
    }
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            mat[i][j]=0;
            for (l = 0; l < dim; ++l) {
                mat[i][j] += matrix1[i][l] * matrix2[l][j];
          }
        }
    }
    return mat;
}
/*
 * calculate the number of col in fileName
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

/*
 * Calculate and output the Weighted Adjacency Matrix
 */
void printWAM(double** matrix, int dim, int num) {
    double **wMatrix = getWeightedMatrix(matrix, dim, num);
    printMat(wMatrix, num);
    freeMemory(wMatrix, num);
}
/*
 * Calculate and output the Diagonal Degree Matrix
 */
void printDDG(double** matrix, int dim, int num) {
    double **dMatrix = getDiagonalDegreeMatrix(matrix, dim,num);
    printMat(dMatrix, num);
    freeMemory(dMatrix, num);
}

void printLNORM(double** matrix, int dim, int num) {
    double **lMatrix = getLaplacianMatrix(matrix, dim, num);
    printMat(lMatrix, num);
    freeMemory(lMatrix, num);
}
/*
 * Outputs must be formatted to 4 decimal places.
 */
void printMat(double** matrix, int dim){
    int i, j;
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            printf("%.4f", matrix[i][j]);
        }
        printf("\n");
    }
}

void printJACOBI(double** matrix, int dim) {
    double* eigenvalues;
    double** V;
    double** result;

    if(!isSymmetric(matrix,dim)){
        printf("Invalid Input!");
        exit(0);
    }

    eigenvalues = (double *) malloc(dim* sizeof (double));
    if(!eigenvalues){
        printf("Invalid Input!");
        exit(0);
    }

    V = jacobiAlgorithm(matrix,dim,eigenvalues);
    result = concatenation(V,eigenvalues,dim);
    freeMemory(V,dim);
    free(eigenvalues);

    printMat(result,dim);
    freeMemory(result,dim);
}

/*
 *  Jacobi algorithm
 */
double** jacobiAlgorithm(double** A, int dim, double* eigenvalues){
    double** pMatrix;
    double** ptMatrix;
    double** V;
    double** tmpV;
    double** Aprime;
    double** tmpA;
    int i;
    int itr = 0;

    V = getUnitMatrix(dim);

    /*
     * Repeat (a),(b) : until A' is diagonal matrix
     *                  until maximum number of rotations = 100
     *                  until the convergence < pow(10,-5)
     */
    while (itr < MAX_ITER){
        /*
         * (a) building a rotation matrix P,PT
         */
        pMatrix = getRotationMat(A,dim);
        ptMatrix = transpose(pMatrix,dim);
        /*
         * V = P1 * P2 * . . .
         */
        tmpV = multiplyMatrices(V,pMatrix,dim);
        freeMemory(V,dim);
        V = tmpV;
        /*
         * (b) Transform the matrix A to: A' = P^T * A * P
         */
        tmpA = multiplyMatrices(A,pMatrix,dim);
        Aprime = multiplyMatrices(ptMatrix,tmpA,dim);

        if(convergence(A,Aprime,dim)){
            freeMemory(A,dim);
            A = Aprime;
            freeMemory(pMatrix,dim);
            freeMemory(ptMatrix,dim);
            break;
        }

        A = Aprime;
        itr++;

        freeMemory(tmpA,dim);
        freeMemory(pMatrix,dim);
        freeMemory(ptMatrix,dim);
    }

    for (i = 0; i < dim; i++) {
        eigenvalues[i] = A[i][i];
    }

    freeMemory(tmpA,dim);
    return V;
}

/*
 * check if off(A)^2 - off(A')^2 <= EPSILON
 */
bool convergence(double** matrix1,double** matrix2,int dim){
    double num = offMatrix(matrix1,dim)- offMatrix(matrix2,dim);
    if(num <= EPSILON){
        return true;
    }
    return false;
}
/*
 * Let off(A)^2 be the sum of squares of all off-diagonal elements of A respectively.
 */
double offMatrix(double** matrix, int dim) {
    double value = 0;
    int i, j;
    for (i = 0; i < dim; i++) {
        for (j = i; j < dim; j++) {
            if (i != j) {
                value += pow(matrix[i][j], 2.0) + pow(matrix[j][i],2);
            }
        }
    }
    return value;
}

/*
 * Rotation Matrix P
 */
double** getRotationMat(double** A, int dim){
    double** pMatrix;
    int i, j, sign_theta;
    double maxValue;
    int maxROW, maxCOL;
    double theta;
    double t, c, s;

    maxValue = A[0][1];
    maxROW = 0;
    maxCOL = 1;
    pMatrix = getUnitMatrix(dim);

    /*
     * Pivot :  The Aij is the off-diagonal element with the largest absolute value.
     *          Aij = maxValue , i = maxROW, j = maxCOL.
     */
    for (i = 0; i < dim; i++) {
        for ( j = i; j < dim; j++) {
            if (fabs(maxValue) < fabs(A[i][j])) {
                maxValue = A[i][j];
                maxROW = i;
                maxCOL = j;
            }
        }
    }

    /*
     * Obtain c and t.
     */
    theta = ((A[maxCOL][maxCOL] - A[maxROW][maxROW])
             / (2*A[maxROW][maxCOL]));
    sign_theta = getSign(theta);
    t = (sign_theta / (fabs(theta) + sqrt(pow(theta, 2) + 1)));
    c = (1 / sqrt(pow(t, 2) + 1));
    s = t * c;


    pMatrix[maxROW][maxROW] = c;
    pMatrix[maxROW][maxCOL] = s;
    pMatrix[maxCOL][maxCOL] = c;
    pMatrix[maxCOL][maxROW] = -s;

    return pMatrix;
}

double** concatenation(double **V, double* eigenvalues, int dim){
    int i,j;
    double **result;
    result = createMat(dim+1, dim);
    if(!result){
        printf("Invalid Input!");
        exit(0);
    }
    for(i=0;i<dim+1;i++){
        for(j=0;j<dim;j++){
            if(i==0){
                result[i][j] = eigenvalues[j];
            }
            else result[i][j] = V[i-1][j];
        }
    }
    return result;
}
