#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "assert.h"
#include <string.h>
#include "spkmeans.h"
#include <stdbool.h>


/*
 * argc := number of inputs.
 * argv := ["...", k, goal, input_filename].
 * argv[0] is the name of the program.
 */
int main(int argc, char* argv[]) {
    FILE *input_file;
    int k, i, n, d, j;
    char* goal;
    Goal g;
    double **inputMatrix;
    int validGoal =0;

    assert(argc > 0);
    if (argc != 4) {
        printf("Invalid Input!");
        return 0;
    }

    /*
     * read the input
     */ 
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
     * check if d==n?
     */
    n = calculateRows(input_file);  /* columns */
    d = calculateCol(input_file);   /* rows */

    /*
     * build a matrix from input file
     */
    inputMatrix = createMat(n, d);
    if (inputMatrix == NULL) {
        printf("Invalid Input!");
        return 0;
    }
    fillMat(input_file, inputMatrix);


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
    opertaion(inputMatrix, k, n, g);
    fclose(input_file);
    freeMemory(inputMatrix, n);
    return 0;
}


void operation(double **matrix, int k, int dim, Goal g) {
    if (g == SPK) {
        double **clusters = createMat(dim, dim);
        spkMeans(EPSILON, matrix, clusters);
    }

    if (g == WAM) {
        printWAM(matrix, dim);
    }
    if (g == DDG) {
        printDDG(matrix, dim);
    }
    if (g == LNORM) {
        printLNORM(matrix, dim);
    }
    if (g == JACOBI) {
        printJACOBI(matrix, dim);
    }
}


double** spkMeans(double epsilon, double** inputMat, double** clusters) {
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

void algorithm(double** clusters, double** inputMat, double** GroupOfClusters, double epsilon){
    int numOfItr=0;
    int x_i;
    int m_i;
    int index;
    while (numOfItr < MAX_ITER){
        if(normMat(clusters, epsilon)==1){
            break;
        }
        resetMat(k,n,GroupOfClusters);
        /*
         * for xi, 0 < i ≤ N:
         *  Assign xi to the closest cluster Sj : argmin_Sj(xi − µj )^2 , ∀j 1 ≤ j ≤ K
         *  0 < index ≤ K
         */
        for(x_i=0;x_i<n;x_i++){
            index = minIndex(clusters, inputMat[x_i]);
            GroupOfClusters[index][x_i] = 10;
        }
        /*
         * for µk, 0 < i ≤ K:
         *  we want to change clusters[i] , 0< i <= k
         */
        for(m_i=0;m_i<k;m_i++){
            update(clusters[m_i], GroupOfClusters[m_i], inputMat);
        }
        numOfItr++;
    }

}


int normMat(double** matrix, double epsilon){
    int i;
    for(i=0;i<k;i++){
        if(calculateNorm(matrix[i]) > epsilon){
            return 0;
        }
    }
    return 1;
}


double calculateNorm(double* vector){
    int i;
    double sum=0.0;
    double value;
    for(i=0;i<d;i++){
        sum = sum + pow(vector[i],2.0);
    }
    value = sqrt(sum);
    return value;
}


void resetMat(int row,int col,double** mat){
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            mat[i][j] = 0.0;
        }
    }
}


void sumVectors(double* vector1,double* vector2){
    int i;
    for(i=0;i<d;i++){
        vector1[i] = vector1[i] + vector2[i];
    }
}


void avgVectors(double* vector1,int cnt){
    int i;
    if(cnt == 0){
        return;
    }
    for(i=0;i<d;i++){
        vector1[i] = vector1[i]/cnt;
    }
}


void update(double* toChage,double* GroupOfClusters,double** inputMat){
    int i;
    int cnt=0;
    /*
     * fill tochange with 0.0
     */
    for(i=0;i<d;i++){
        toChage[i]= 0.0;
    }
    for(i=0;i<n;i++){
        if(GroupOfClusters[i] > 5){
            cnt++;
            sumVectors(toChage,inputMat[i]);
        }
    }
    avgVectors(toChage,cnt);
}


int minIndex(double** clusters, double* victor){
    int minIndex=0;
    double minDistance=DBL_MAX;
    double tempMinDistance;
    int i;

    for(i=0;i<k;i++){
        tempMinDistance = distance(victor, clusters[i]);
        if(tempMinDistance<minDistance) {
            minIndex = i;
            minDistance = tempMinDistance;
        }
    }
    return minIndex;
}


double distance(double* vector1 , double* vector2){
    int i;
    double sum=0.0;
    for(i=0; i < d; i++){
        /*
         * argmin_Sj(xi − µj)^2
         */
        sum = sum + pow((vector1[i] - vector2[i]), 2.0);
    }
    return sum;
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
double** getWeightedMatrix(double** matrix, int dim){
    int i, j;
    double** wMatrix;
    wMatrix = createMat(dim, dim);
    if(!wMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    for (i = 0; i < dim; i++) {
        for (j = (i+1); j < dim; j++) {
            wMatrix[i][j] = calculateWeight(matrix,i , j, dim);
            wMatrix[j][i] = wMatrix[i][j];
        }
    }
    return wMatrix;
}

/*
 * if i=j dij = sumOf(wiz) for all z=1..n, otherwise dij=0
 */
double** getDiagonalMatrix(double** matrix, int dim){
    int i, j;
    double sum;
    double** dMatrix;
    dMatrix = createMat(dim, dim);
    if(!dMatrix){
        printf("An Error Has Occurred");
        exit(0);
    }
    for (i = 0; i < dim; ++i){
        sum = 0;
        for (j = 0; j < dim; j++) {
            sum += matrix[i][j];
        }
        dMatrix[i][i] = sum;
    }
    return dMatrix;
}

/*
 * Form the diagonal degree matrix D ∈ R^(n×n).
 * the diagonal equals to the sum of the i-th row of W.
 */
double** getDiagonalDegreeMatrix(double** matrix, int dim){
    int i;
    double** dMatrix;
    double** wMatrix;
    wMatrix = getWeightedMatrix(matrix, dim);
    dMatrix = getDiagonalMatrix(wMatrix, dim);
    for (i = 0; i < dim; ++i) {
        dMatrix[i][i] = pow(sqrt(dMatrix[i][i]), -0.5);
    }
    freeMemory(wMatrix,dim);
    return dMatrix;
}

/*
 * Form The Normalized Graph Laplacian matrix Lnorm ∈ R^(n×n).
 * Lnorm = I − (D^(-1/2) * W * D^(-1/2))
 *       = ( (Diagonal Degree Matrix * Weighted Matrix) * Diagonal Degree Matrix )
 */
double** getLaplacianMatrix(double** matrix, int dim){
    int i, j;
    double** wMatrix;
    double** dMatrix;
    double** lnormMatrix;
    double** temp;
    wMatrix = getWeightedMatrix(matrix,dim);
    if(!wMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    dMatrix = getDiagonalDegreeMatrix(matrix,dim);
    if(!dMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    /*
     * temp = D^(-1/2) * W
     */
    temp = multiplyMatrices(dMatrix, wMatrix, dim);
    if(!temp){
        printf("An Error Has Occured");
        exit(0);
    }
    /*
     * lnormMatrix = (D^(-1/2) * W)*D^(-1/2)
     */
    lnormMatrix = multiplyMatrices(temp, dMatrix, dim);
    if(!lnormMatrix){
        printf("An Error Has Occured");
        exit(0);
    }
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (i == j) {
                lnormMatrix[i][j] = 1 - lnormMatrix[i][j];
            } else {
                lnormMatrix[i][j] = ((-1)*(lnormMatrix[i][j]));
            }
        }
    }
    freeMemory(temp, dim);
    freeMemory(wMatrix,dim);
    freeMemory(dMatrix,dim);
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
 * Let off(A)^2 be the sum of squares of all off-diagonal elements of A respectively.
 */
double offOfMat(double** matrix, int dim) {
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
        for ( j = i + 1; j < dim; j++) {
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
    double **pMatrix = getRotationMat(matrix, dim);      //TODO: Check if we are in even iteration?
    double **vMatrix = getRotationMat(matrix, dim);
    // building matrix A'
    double **aTagTemp = multiplyMatrices(transpose(pMatrix, dim), matrix, dim);
    double **aTagMat = multiplyMatrices(aTagTemp, pMatrix, dim);

    offset = offOfMat(matrix, dim) - offOfMat(aTagMat, dim);

    // repeat until A' is diagonal matrix
    while ((isDiagonalMatrix(aTagMat, dim) == 0) || (iter < 100 || EPSILON < offset )) {
        pMatrix = getRotationMat(aTagMat, dim);
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
 *  matrix full of zeroz
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
/*
 *  multiplying matrices dimXdim
 *  the number of columns in the first matrix must be equal
 *  to the number of rows in the second matrix.
 */
double** multiplyMatrices(double** matrix1, double** matrix2, int dim) {
    int j, i, k;
    double** mat;
    mat = createMat(dim, dim);
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            mat[i][j]=0;
            for (k = 0; k < dim; ++k) {
                mat[i][j] += matrix1[i][k] * matrix2[k][j];
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


void printWAM(double** matrix, int dim) {
    double **wMatrix = getWeightedMatrix(matrix, dim);
    printMat(wMatrix, dim);
    freeMemory(wMatrix, dim);
}

void printDDG(double** matrix, int dim) {
    double **wMatrix = getWeightedMatrix(matrix, dim);
    double **dMatrix = getDiagonalD(matrix, dim);
    printMat(dMatrix, dim);
    freeMemory(wMatrix, dim);
    freeMemory(dMatrix, dim);
}

void printLNORM(double** matrix, int dim) {
    double **lMatrix = getLaplacianMatrix(matrix, dim);
    printMat(lMatrix, dim);
    freeMemory(lMatrix, dim);
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
void printJACOBI(double** matrix, int dim) {}