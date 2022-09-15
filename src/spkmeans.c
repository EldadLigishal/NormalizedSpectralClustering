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
     *  close file and freeing memory
     */
    fclose(ifp);
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

int getSign(double number){
    if (number >= 0) {
        return 1;
    }
    return -1;
}

double** getUnitMatrix(int dim){
    int i;
    double ** matrix = createMat(dim, dim);
    if(!matrix){
        printf("Invalid Input!");
        exit(0);
    }
    for ( i = 0; i < dim; i++) {
        matrix[i][i] = 1;
    }
    return matrix;
}

double** createMat(int col, int row){
    int i, j;
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
    for (i = 0; i < col; i++){
        for (j = 0; j < row; j++){
            matrix[i][j] = 0;
        }
    }
    return matrix;
}

void printMatJacobi(double** matrix, int dim, int num){
    int i;
    for(i=0;i<dim;i++){
        if(matrix[0][i]<0 && matrix[0][i] > - 0.00001){ matrix[0][i]=0 ; }
    }
    printMat(matrix,num+1,dim);
    freeMemory(matrix,num+1);
}

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
    while(token != NULL){
        token = strtok(NULL, breaks);
        ++cnt;
    }
    /*
     * close file
     */
    fclose(ifp);
    return cnt;

}

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
        while (token != NULL) {
            inputMat[row][col] = strtod(token, &useless);
            token = strtok(NULL,breaks);
            col++;
        }
        col=0;
        row++;
    }
    fclose(ifp);
}

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

void printWAM(double** matrix, int dim, int num){
    double **wMatrix;
    wMatrix = getWeightedMatrix(matrix, dim, num);
    printMat(wMatrix, num, num);
    freeMemory(wMatrix, num);
}

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

void printMat(double** matrix, int rows, int cols){
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            printf("%.4f",matrix[i][j]);
            if(j!=cols-1) printf(",");
        }
        printf("\n");
    }
}

void printJACOBI(double** matrix, int dim, int num){
    double* eigenvalues;
    double** V;
    double** result;
    int i;

    eigenvalues = (double *) malloc(num* sizeof (double));
    if(!eigenvalues){
        printf("Invalid Input!");
        exit(0);
    }
    for (i = 0; i < num; i++) {
        eigenvalues[i] = 0.0;
    }

    V = jacobiAlgorithm(matrix,num,eigenvalues);
    result = concatenation(V,eigenvalues,num);

    printMatJacobi(result,dim,num);
    freeMemory(V,num);
    free(eigenvalues);
}

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

bool convergence(double** matrix1,double** matrix2,int dim){
    double num = offMatrix(matrix1,dim)- offMatrix(matrix2,dim);
    if(num <= EPSILON){
        return true;
    }
    return false;
}

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

void RotationMat(double** A,double** pMatrix,double** ptMatrix,int dim){
    int i, j, sign_theta;
    double maxValue;
    int maxROW, maxCOL;
    double theta;
    double t, c, s;
    maxValue = A[0][1];
    maxROW = 0;
    maxCOL = 1;
    for (i = 0; i < dim; i++){
        pMatrix[i][i] = 1.0;
        ptMatrix[i][i] = 1.0;
    }
    /*
     * Pivot :  The Aij is the off-diagonal element with the largest absolute value.
     *          Aij = maxValue , i = maxROW, j = maxCOL.
     */
    for (i = 0; i < dim; i++){
        for ( j = i + 1; j < dim; j++){
            if (fabs(maxValue) < fabs(A[i][j]) && (i != j)){
                maxValue = A[i][j];
                maxROW = i;
                maxCOL = j;
            }
        }
    }
    /*
     * Obtain c and t.
     */
    if (A[maxROW][maxCOL] != 0) {
        theta = ((A[maxCOL][maxCOL] - A[maxROW][maxROW])
                 / (2*A[maxROW][maxCOL]));
    } else {
        theta = 0.0;
    }
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
            } else {
                result[i][j] = V[i-1][j];
            }
        }
    }
    return result;
}

int getEigengapHeuristic(double* array,int len){
    int maxIndex;
    int i;
    double deltaMax;
    deltaMax = fabs(array[0] - array[1]);
    maxIndex = 0;
    /*
     * δ_max >= δi (= |λi − λi+1|)
     */
    for (i = 0; i < floor(len / 2); i++){
        if (fabs(array[i] - array[i + 1]) > deltaMax) {
            maxIndex = i;
            deltaMax = fabs(array[i] - array[i + 1]);
        }
    }
    return maxIndex+1;
}

double** getTMatrix(double** matrix,int dim, int num, int k){
    int i,j,*index;
    double sq,sum;
    double **lMatrix,**vMatrix,**tmp1,**tmp2,**uMatrix,**tMatrix,*eigenvalues;
    tmp1=createMat(num, num);
    tmp2=createMat(num, num);
    index = (int*) malloc(num * (sizeof(int)));
    if(!index){
        printf("An Error Has Occurred\n");
        exit(0);
    }
    for(i=0;i<num;i++){
        index[i]=i;
    }
    /*
     * 1: Form the weighted adjacency matrix W from X(X:=matrix)
     * 2: Compute the normalized graph Laplacian Lnorm(Lnorm:=lMatrix)
     */
    lMatrix = getLaplacianMatrix(matrix, dim, num);
    /*
     * 3: Determine k and obtain the largest k eigenvectors u1,...,uk of Lnorm,
     * eigenvalues must be ordered decreasingly.
     * Determining k would be based on the eigengap heuristic or given as an input.
     */
    eigenvalues = (double*) malloc(num * sizeof(double));
    if (!eigenvalues){
        printf("An Error Has Occurred\n");
        exit(0);
    }
    for (i =0 ; i < num ; i++) {
        eigenvalues[i] = 0.0;
    }
    vMatrix = jacobiAlgorithm(lMatrix,num,eigenvalues);
    sortingEigenValues(eigenvalues, num, index);
    if(k==0){
        K = getEigengapHeuristic(eigenvalues,num);
        k = K;
    } else {
        K = k;
    }
    if(k>num){
        printf("An Error Has Occurred\n");
        exit(0);
    }
    /*
     * 4: Let U ∈ Rn×k be the matrix containing the vectors u1,...,uk as columns
     */
    updateTmp1(tmp1,vMatrix,num);
    updateTmp2(tmp2,tmp1,index,num);
    uMatrix = createMat(num, k);
    tMatrix = createMat(num, k);
    updateUmatrix(uMatrix,tmp2,num,k);
    /*
     * 5: Form the matrix T ∈ Rn×k from U by renormalizing each of U’s rows to have unit length,
     *     that is set ...
     */
    for(i=0;i<num;i++){
        sum=0.0;
        for(j=0;j<k;j++){
            sum+=pow(uMatrix[i][j],2);
        }
        sq= sqrt(sum);
        for(j=0;j<k;j++){
            if(sq!=0.0){
                tMatrix[i][j] = (uMatrix[i][j]) / sq;
            }else{
                tMatrix[i][j]=(uMatrix[i][j]);
            }
        }
    }
    freeMemory(vMatrix, num);
    freeMemory(tmp1, num);
    freeMemory(tmp2, num);
    freeMemory(uMatrix, num);
    free(eigenvalues);
    free(index);
    return tMatrix;
}
void updateTmp1(double ** tmp1,double **matrix,int dim){
    int i,j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            tmp1[i][j]=matrix[j][i];
        }
    }
}
void updateTmp2(double ** tmp2,double **tmp1,int* index,int dim){
    int i,j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            tmp2[i][j]=tmp1[index[i]][j];
        }
    }
}
void updateUmatrix(double **uMatrix,double **tmp2,int num,int k){
    int i,j;
    for(i=0;i<num;i++) {
        for(j=0;j<k;j++){
            uMatrix[i][j]=tmp2[j][i];
        }
    }
}

void sortingEigenValues(double *arr_double, int len,int* arr_int) {
    int i, j, temp_int;
    double temp_double;
    for (i = 0; i < len; i++) {
        for (j = i + 1; j < len; j++) {
            if (arr_double[i] < arr_double[j]) {
                temp_double = arr_double[i];
                temp_int = arr_int[i];

                arr_double[i] = arr_double[j];
                arr_int[i] = arr_int[j];

                arr_double[j] = temp_double;
                arr_int[j] = temp_int;
            }
        }
    }
}

void kmeans(int k , int maxItr, int d, int num, double **inputMatrix, double **init_centroids){
    int vector_dim, i, j, centroid_index, iteration, all_vectors_num;
    double **sum_mat, **old_centroids;
    int *clusters_size;
    all_vectors_num=num;
    vector_dim = d;
    iteration = 0;
    centroids = createMat(k,d);
    if (!centroids){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i=0;i<k;i++){
        for (j=0;j<d;j++) {
            centroids[i][j] = init_centroids[i][j];
        }
    }
    freeMemory(init_centroids,k);
    sum_mat = createMat(k,vector_dim);
    /*
     * initializing clusters sizes to zero
     */
    clusters_size = (int *) malloc(k * sizeof(int));
    if (!clusters_size) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    old_centroids = createMat(k,vector_dim);
    if (!old_centroids) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    /*
     * ALGORITHM LOOP
     */
    while (iteration < maxItr) {
        /* reset cluster sizes and sum matrix */
        for (i = 0; i < k; i++) {
            clusters_size[i] = 0;
            for (j = 0; j < vector_dim; j++) {
                sum_mat[i][j] = 0;
            }
        }
        for (j = 0; j < all_vectors_num; j++) {
            centroid_index = calculateNorm(centroids, inputMatrix[j], k, vector_dim);/*calc norm returns index of the current closest centroid to the vector*/
            clusters_size[centroid_index] += 1;
            for (i = 0; i < vector_dim; i++) {
                sum_mat[centroid_index][i] += inputMatrix[j][i];
            }
        }
        for (i = 0; i < k; i++) {
            for (j = 0; j < vector_dim; j++) {
                old_centroids[i][j] = centroids[i][j];
            }
        }
        for (i = 0; i < k; i++) {
            free(centroids[i]);
            centroids[i] = divide(sum_mat[i], clusters_size[i], vector_dim);
        }
        iteration += 1;
    }
    freeMemory(inputMatrix, all_vectors_num);
    freeMemory(old_centroids, k);
    freeMemory(sum_mat,k);
    free(clusters_size);
}

double* divide(double* vector,int num,int len){
    double *result;
    double number;
    int i;
    result = (double*) malloc(len * sizeof(double));
    if(!result){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < len; i++) {
        number = vector[i]/num;
        result[i] = number;
    }
    return result;
}

int calculateNorm(double** matrix, double* vector, int k, int dim){
    double normValue;
    double tmp;
    int i;
    int value = 0;
    normValue = calculateNormDistance(matrix[0], vector, dim);
    for(i=1;i<k;i++){
        tmp = calculateNormDistance(matrix[i], vector, dim);
        if(tmp < normValue){
            value = i;
            normValue = tmp;
        }
    }
    return value;
}
double calculateNormDistance(double* vector1, double* vector2, int len){
    double* arr;
    int i;
    double sum = 0;
    double value;
    arr = (double*) malloc(len * sizeof (double));
    if(!arr){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i=0; i < len; i++){
        value = vector1[i]-vector2[i];
        arr[i]= value;
    }
    for(i=0; i < len; i++){
        value = sum + arr[i]*arr[i];
        sum = value;
    }
    free(arr);
    return sum;
}