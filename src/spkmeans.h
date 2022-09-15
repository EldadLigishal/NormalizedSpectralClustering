#ifndef SPKMEANS_LIBRARY_H
#define SPKMEANS_LIBRARY_H

#include <stdbool.h>
#define EPSILON 0.00001
#define LINESIZE 1000
#define MAX_ITER 100

/**
 * global variable named centroids .
 * this variable will help us in calculating the centroids when goal=spk
 */
double **centroids;
/**
 * global variable named K = the number of clusters .
 */
int K;


/**
 * Goal enums in C file.
 */
typedef enum {WAM, DDG, LNORM, JACOBI
} Goal;

/**
 *
 * @param argc number of inputs.
 * @param argv ["...", goal,filename].
 * @return
 */
int main(int argc, char* argv[]);
/**
 *
 * @param matrix matrix of dimension dimXnum
 * @param dim is the number of column of an input file.
 * @param num is the number of line of an input file , = number of vectors.
 * @param g is the enum which include the type if algorithm we have to run.
 * determine which function to run depending on the goal
 */
void operation(double **matrix, int dim, int num, Goal g);
/**
 *
 * @param matrix matrix of dimension rowsXcols.
 * @param rows is the number of line of given matrix , = number of vectors.
 * @param cols is the number of column of given matrix.
 * print a 2-dimension matrix
 * Outputs must be formatted to 4 decimal places.
 */
void printMat(double** matrix, int rows, int cols);
/**
 *
 * @param matrix matrix of dimension dimXnum.
 * @param dim is the number of column of given matrix.
 * @param num is the number of line of given matrix , = number of vectors.
 * calculates the Weighted Adjacency Matrix and print it.
 */
void printWAM(double** matrix, int dim, int num);
/**
 *
 * @param matrix matrix of dimension dimXnum
 * @param dim is the number of column of given matrix.
 * @param num is the number of line of given matrix , = number of vectors.
 * calculates the Diagonal Degree Matrix and print it.
 */
void printDDG(double** matrix, int dim, int num);
/**
 *
 * @param matrix matrix of dimension dimXnum
 * @param dim is the number of column of given matrix.
 * @param num the number of line of given matrix , = number of vectors.
 * calculates the Laplacian Matrix and print it.
 */
void printLNORM(double** matrix, int dim, int num);
/**
 *
 * @param matrix matrix of dimension dimXnum
 * @param dim is the number of column of given matrix.
 * @param num is the number of line of given matrix , = number of vectors.
 * calculates the Jacobi Matrix and print it.
 */
void printJACOBI(double** matrix, int dim, int num);
/**
 *
 * @param matrix matrix with dimension dimXnum
 * @param dim is the number of column of an input matrix.
 * @param num is the number of line of an input matrix , = number of vectors.
 * prints Jacobi matrix
 */
void printMatJacobi(double** matrix, int dim, int num);

/**
 *
 * @param matrix square matrix of dimension dimXdim
 * @param index1 a number that represent index
 * @param index2 a number that represent index
 * @param dim a number that represent the dimension of input matrix
 * @return the value of exp(−||xi − xj||/2)
 * calculates the value of exp(−(||xi − xj||/2)), and returns it.
 */
double calculateWeight(double** matrix, int index1, int index2, int dim);
/**
 *
 * @param matrix matrix of dimension dimXnum
 * @param dim is the number of column of an input file.
 * @param num is the number of line of an input file , = number of vectors.
 * @return weighted adjacency matrix
 * Form the weighted adjacency matrix W ∈ R^(n×n) , and returns it.
 * The weights are symmetric (wij = wji) and non-negative (wij ≥ 0).
 * wii = 0 for all i’s
 * the rest of the values are set to: wij = exp(−||xi − xj||/2)
 */
double** getWeightedMatrix(double** matrix, int dim, int num);
 /**
  *
  * @param matrix matrix of dimension dimXnum
  * @param dim is the number of column of an input matrix.
  * @param num is the number of line of an input matrix , = number of vectors.
  * @return the Diagonal Degree matrix
  * Form the diagonal degree matrix D ∈ R^(n×n) , and returns it.
  * the diagonal equals to the sum of the i-th row of W.
  * if the index i=j then dij=sumOf(wiz) for all z=1..n, otherwise dij=0.
 */
double** getDiagonalDegreeMatrix(double** matrix,int dim, int num);
 /**
  *
  * @param matrix matrix of dimension dimXnum
  * @param dim is the number of column of an input matrix.
  * @param num is the number of line of an input  matrix, = number of vectors.
  * @return the Normalized Graph Laplacian matrix
  * Form The Normalized Graph Laplacian matrix Lnorm ∈ R^(n×n) , and returns it.
  * Lnorm = I − (D^(-1/2) * W * D^(-1/2))
  *       = ( (Diagonal Degree Matrix * Weighted Matrix) * Diagonal Degree Matrix )
 */
double** getLaplacianMatrix(double** matrix, int dim, int num);
/**
 *
 * @param A a symmetric matrix
 * @param n is the number of line of input matrix , = number of vectors.
 * @param eigenvalues array of eigenvalues numbers
 * @return the Jacobi matrix
 * Calculate and output the eigenvalues and eigenvectors
 * performs the jacobi algorithm that described in the project instructions
*/
double** jacobiAlgorithm(double** A, int n, double* eigenvalues);
/**
 *
 * @param dim the dimension of the matrix
 * @return the Identity Matrix
 * this function receives a number that represent the dimension of the Identity matrix,
 * and returns the identity matrix.
*/
double** getUnitMatrix(int dim);
/**
 *
 * @param A a symmetric matrix
 * @param pMatrix  rotation matrix P
 * @param ptMatrix the transpose of rotation matrix P
 * @param dim is the number of column of A.
 * builds a rotation matrix P,PT.
 */
void RotationMat(double** A,double** pMatrix,double** ptMatrix,int dim);
/**
 *
 * @param number the number that we want to check its sign
 * @return 1 if x>=0, else 0
 */
int getSign(double number);
/**
 *
 * @param matrix1  matrix with dimension = dimXdim
 * @param matrix2  matrix with dimension = dimXdim
 * @param dim the dimension of matrix1 and matrix2
 * @return  matrix1 * matrix2
 * performs matrix multiplication.
 * the number of columns in the first matrix must be equal
 * to the number of rows in the second matrix.
 */
double** multiplyMatrices(double** matrix1, double** matrix2, int dim);
/**
 * @param matrix matrix with dimension = dimXdim
 * @param dim the dimension of matrix
 * @return
 * Let off(A)^2 be the sum of squares of all off-diagonal elements of A respectively.
 */
double offMatrix(double** matrix, int dim);
/**
 *
 * @param matrix1 matrix with dimension = dimXdim
 * @param matrix2 matrix with dimension = dimXdim
 * @param dim the dimension of matrix1 and matrix2
 * @return
 * check if off(A)^2 - off(A')^2 <= EPSILON
 * computes the convergence condition to stop the jacobi algorithm
 */
bool convergence(double** matrix1,double** matrix2,int dim);
/**
 *
 * @param V the Jacobi matrix
 * @param eigenvalues array of eigenvalues numbers
 * @param dim is the number of column of an input file.
 * @return
 */
double** concatenation(double **V, const double *eigenvalues, int dim);
/**
 * @param matrix matrix with number of columns = len.
 * @param len number of the columns that are in the given matrix.
 * free memory function.
*/
void freeMemory(double** matrix ,int len);
/**
 * @param col number of column in matrix
 * @param row number of rows in matrix
 * @return matrix with dimension col x row
 * creates and returns a new matrix with dimension col x row
*/
double** createMat(int col, int row);
/**
 * @param fileName *.txt file that contains data-points separated by commas.
 * @return the number of column that are in fileName.
 * calculate the number of column in fileName and returns it.
 */
int calculateCol(char* fileName);
/**
 * @param fileName *.txt file that contains data-points separated by commas.
 * @return the number of rows that are in fileName.
 * calculate the number of rows in fileName and returns it.
 */
int calculateRows(char* fileName);
/**
 * @param fileName *.txt file that contains data-points separated by commas.
 * @param inputMat matrix with dimension = (number of rows in given file)X(number of col in given file).
 * fill 2-dimensional arrays.
 * reads the input points from the given file,
 * and updates the given matrix to contain the points.
 */
void fillMat(char* fileName,double** inputMat);

/**
 *
 * @param row is the number of line of given matrix , = number of vectors.
 * @param col is the number of column of given matrix.
 * @param mat matrix with dimension = rowXcol.
 * reset all matrix values to 0
 */
void resetMat(int row,int col,double** mat);

/**
 *
 * @param array array of eigenvalues numbers.
 * @param len the length of given array.
 * @return the number of clusters k.
 * The Eigengap Heuristic
 * In order to determine the number of clusters k,
 * we will use eigengap heuristic.
 */
int getEigengapHeuristic(double* array,int len);

/**
 * @param matrix matrix with dimension = dimXnum
 * @param dim is the number of column of given matrix.
 * @param num is the number of line of given matrix , = number of vectors.
 * @param k the number of clusters k.
 * @return T matrix
 * performs The Normalized Spectral Clustering Algorithm
 * as desctibed in the first page of the project instructions and returns
 * the T matrix
 */
double** getTMatrix(double** matrix,int dim, int num, int k);
/**
 *
 * @param tmp1 matrix with dimension = dimXnum
 * @param matrix matrix with dimension = dimXnum
 * @param dim  is the number of column of given matrix.
 * update the tmp1 matrix in getTMatrix function
 */
void updateTmp1(double ** tmp1,double **matrix,int dim);
/**
 *
 * @param tmp2 matrix with dimension = dimXdim
 * @param tmp1 matrix with dimension = dimXdim
 * @param index array of index
 * @param dim is the number of column of given matrix.
 * update the tmp2 matrix in getTMatrix function.
 */
void updateTmp2(double ** tmp2,double **tmp1,int* index,int dim);
/**
 *
 * @param uMatrix matrix with dimension = dimXdim
 * @param tmp2 matrix with dimension = numXnum
 * @param num is the number of rows of given matrix.
 * @param k the number of clusters k
 * update the uMatrix matrix in getTMatrix function.
 */
void updateUmatrix(double **uMatrix,double **tmp2,int num,int k);
/**
 *
 * @param k the number of clusters k
 * @param maxItr = 300
 * @param d is the number of column of given matrix.
 * @param num is the number of rows of given matrix.
 * @param inputMatrix matrix with dimension = dXnum
 * @param init_centroids matrix with dimension = dXnum
 * This function updates the clusters (performs steps 5 and 6 f`rom hw2)
*/
void kmeans(int k , int maxItr, int d, int num, double **inputMatrix, double **init_centroids);
/**
 *
 * @param vector array of numbers
 * @param num divides every number in vector with this num
 * @param len the length of given vector
 * @return a new vector = array of numbers
 * divide every number in the given vector with given num.
 */
double* divide(double* vector,int num,int len);
/**
 *
 * @param vector1 array of numbers
 * @param vector2 array of numbers
 * @param len the length of given vector
 * @return distance between tow vectors
 * calculate the distance between two vectors , and returns it.
*/
double calculateNormDistance(double* vector1, double* vector2, int len);
/**
 *
 * @param matrix matrix with dimension = dimXdim
 * @param vector array of numbers.
 * @param k the number of clusters k
 * @param dim is the number of column of given matrix.
 * @return
 * This function calculates and returns the euclidean distance between two vectors.
 */
int calculateNorm(double** matrix, double* vector, int k, int dim);
/**
 * @param arr_double array of eigenvalues numbers.
 * @param len the length of given array.
 * @param arr_int array of index
 * sorting the eigenvalues in decreasing order.
 */
void sortingEigenValues(double *arr_double, int len,int* arr_int);
#endif