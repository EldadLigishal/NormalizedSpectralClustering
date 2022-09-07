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
int main(int argc, char* argv[]){
    char* filename;
    FILE* ifp;
    int n, d;
    char* goal;
    Goal g;
    /*
     * input_filename := *.txt file that contains data-points separated by commas.
     */
    double **inputMatrix;
    int validGoal = 0;

    if (argc != 3){
        printf("Invalid Input!");
        return 0;
    }
    /*
     *  open file
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
     *  calculating the columns and rows of the input file
     *  d := number of column of an input file.
     *  n := number of line of an input file , = number of vectors.
     */
    n = calculateRows(filename);
    d = calculateCol(filename);

    /*
     *  build a matrix from input file.
     *  We start by creating an empty matrix, and then we will fill it.
     */
    inputMatrix = createMat(n, d);
    if (inputMatrix == NULL){
        printf("Invalid Input!");
        return 0;
    }
    fillMat(filename, inputMatrix);

    operation(inputMatrix,d,n,g);
    /*
     *  close file
     */
    fclose(ifp);
    /*
     *  freeing memory
     */
    if(g != JACOBI){
        freeMemory(inputMatrix, n);
    }
    return 0;
}

void operation(double **matrix, int dim,int num, Goal g){
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
        printJACOBI(matrix, dim, num);
    }
}
/*
 * calculate the value of exp(−||xi − xj||/2)
 */
double calculateWeight(double** matrix, int index1, int index2, int dim){
    int i;
    double value = 0.0;
    double distance;
    for (i = 0; i < dim; i++){
        distance = fabs(matrix[index1][i] - matrix[index2][i]);
        value += pow(distance, 2.0);
    }
    value = exp(((-1)*( (sqrt(value)) / 2)));
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
        printf("An Error Has Occurred");
        exit(0);
    }
    for (i = 0; i < num; i++) {
        for (j = i; j < num; j++) {
            if(i==j){
                wMatrix[i][j]=0.0;
            }else{
                wMatrix[i][j] = calculateWeight(matrix,i , j, dim);
                wMatrix[j][i] = wMatrix[i][j];
            }
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
        sum = 0.0;
        for (j = 0; j < num; j++){
            sum += wMatrix[i][j];
        }
        dMatrix[i][i] = sum;
    }
    freeMemory(wMatrix,num);
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
    double** I;
    wMatrix = getWeightedMatrix(matrix,dim, num);
    dMatrix = getDiagonalDegreeMatrix(matrix,dim, num);
    I = getUnitMatrix(num);
    /*
     *  calculate D^-0.5
     */
    for(i=0;i<num;i++){
        dMatrix[i][i] = pow(dMatrix[i][i], -0.5);
    }
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
        printf("An Error Has Occurred");
        exit(0);
    }
    for (i = 0; i < num; i++) {
        for (j = 0; j < num; j++) {
            lnormMatrix[i][j] = I[i][j]-lnormMatrix[i][j];
        }
    }
    freeMemory(temp, num);
    freeMemory(I,num);
    freeMemory(wMatrix,num);
    freeMemory(dMatrix,num);
    return lnormMatrix;
}

/*
 * sign(x) = 1 if x >= 0, else 0
 */
int getSign(double number){
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

double** getUnitMatrix(int dim){
    int i, j;
    double ** matrix = createMat(dim, dim);
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
    int i;
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

void printMatJacobi(double** matrix, int dim, int num){
    int i;
    for(i=0;i<dim;i++){
        if(matrix[0][i]<0 && matrix[0][i] > - 0.00005){ matrix[0][i]=0 ; }
    }
    printMat(matrix,num+1,dim);
    freeMemory(matrix,num+1);
}
/*
 *  multiplying matrices dimXdim
 *  the number of columns in the first matrix must be equal
 *  to the number of rows in the second matrix.
 */
double** multiplyMatrices(double** matrix1, double** matrix2, int dim){
    int j, i, l;
    double** mat = NULL;
    mat = createMat(dim, dim);
    assert(mat!=NULL);
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            mat[i][j]=0.0;
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
    if(ifp ==NULL){
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
void printWAM(double** matrix, int dim, int num){
    double **wMatrix;
    wMatrix = getWeightedMatrix(matrix, dim, num);
    printMat(wMatrix, num, num);
    freeMemory(wMatrix, num);
}

/*
 * Calculate and output the Diagonal Degree Matrix
 */
void printDDG(double** matrix, int dim, int num){
    double **dMatrix;
    dMatrix = getDiagonalDegreeMatrix(matrix, dim,num);
    printMat(dMatrix, num, num);
    freeMemory(dMatrix, num);
}

void printLNORM(double** matrix, int dim, int num){
    double **lMatrix;
    lMatrix = getLaplacianMatrix(matrix, dim, num);
    printMat(lMatrix, num, num);
    freeMemory(lMatrix, num);
}

/*
 * Outputs must be formatted to 4 decimal places.
 */
void printMat(double** matrix, int rows, int cols){
    int i, j;
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols-1; j++){
            if(fabs(matrix[i][j])<0.00005){
                printf("%.4f,", fabs(matrix[i][j]));
            }
            else{
                printf("%.4f,", matrix[i][j]);
            }
        }
        if(i==rows-1){
            printf("%.4f", matrix[i][j]);
        }
        else{
            printf("%.4f\n", matrix[i][j]);
        }
    }
}

void printJACOBI(double** matrix, int dim, int num){
    double* eigenvalues;
    double** V;
    double** result;

    eigenvalues = (double *) malloc(num* sizeof (double));
    if(!eigenvalues){
        printf("Invalid Input!");
        exit(0);
    }

    V = jacobiAlgorithm(matrix,num,eigenvalues);
    result = concatenation(V,eigenvalues,num);


    printMatJacobi(result,dim,num);

    freeMemory(V,num);
    free(eigenvalues);
}


/*
 *  Jacobi algorithm
 */
double** jacobiAlgorithm(double** matrix, int n, double* eigenvalues){
    double** pMatrix;
    double** ptMatrix;
    double** V;
    double** tmp;
    double** Aprime;
    double** A;
    double** tmpA;
    int i, itr = 0;

    A = matrix;
    V = getUnitMatrix(n);

    pMatrix = createMat(n,n);
    ptMatrix= createMat(n,n);

    /*
     * Repeat (a),(b) : until A' is diagonal matrix
     *                  until maximum number of rotations = 100
     *                  until the convergence < pow(10,-5)
     */
    while (itr < MAX_ITER){
        resetMat(n,n,pMatrix);
        resetMat(n,n,ptMatrix);
        /*
         * (a) building a rotation matrix P,PT
         */
        RotationMat(A,pMatrix,ptMatrix,n);
        /*
         * V = P1 * P2 * . . .
         */
        tmp = multiplyMatrices(V, pMatrix, n);
        freeMemory(V, n);
        V = tmp;
        /*
         * (b) Transform the matrix A to: A' = P^T * A * P
         */
        tmpA = multiplyMatrices(ptMatrix, A, n);
        Aprime = multiplyMatrices(tmpA, pMatrix, n);
        freeMemory(tmpA,n);

        if(convergence(A, Aprime, n)){
            freeMemory(A, n);
            A = Aprime;
            break;
        }
        freeMemory(A, n);
        A = Aprime;
        itr++;
    }

    for (i = 0; i < n; i++) {
        eigenvalues[i] = A[i][i];
    }
    freeMemory(A,n);
    freeMemory(pMatrix,n);
    freeMemory(ptMatrix,n);

    return V;
}

void resetMat(int row,int col,double** mat){
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            mat[i][j] = 0.0;
        }
    }
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
void RotationMat(double** A,double** pMatrix,double** ptMatrix,int dim){
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
    for (i = 0; i < dim; i++){
        for ( j = i+1; j < dim; j++){
            if (fabs(maxValue) < fabs(A[i][j])){
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

    ptMatrix[maxROW][maxROW] = c;
    ptMatrix[maxROW][maxCOL] = -s;
    ptMatrix[maxCOL][maxCOL] = c;
    ptMatrix[maxCOL][maxROW] = s;

}

double** concatenation(double **V, const double* eigenvalues, int dim){
    int i,j;
    double **result;
    result = createMat(dim+1, dim);
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
