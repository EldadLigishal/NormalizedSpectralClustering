#ifndef SPKMEANS_LIBRARY_H
#define SPKMEANS_LIBRARY_H

#include <stdbool.h>
#define EPSILON pow(10,-5)
#define LINESIZE 1000
#define MAX_ITER 100



typedef enum {WAM, DDG, LNORM, JACOBI
} Goal;

/**
 *
 * @param argc
 * @param argv := ["...", goal,filename].
 * @return
 */
int main(int argc, char* argv[]);
/**
 *
 * @param matrix
 * @param dim is the number of column of an input file.
 * @param num is the number of line of an input file , = number of vectors.
 * @param g is the enum which include the type if algorithm we have to run.
 */
void operation(double **matrix, int dim, int num, Goal g);
/**
 *
 * @param matrix
 * @param dim is the number of column of an input file.
 * @param num is the number of line of an input file , = number of vectors.
 */
void printWAM(double** matrix, int dim, int num);
/**
 *
 * @param matrix
 * @param rows
 * @param cols
 * Outputs must be formatted to 4 decimal places.
 */
void printMat(double** matrix, int rows, int cols);
/**
 *
 * @param matrix
 * @param dim
 * @param num
 */
void printDDG(double** matrix, int dim, int num);
/**
 *
 * @param matrix
 * @param dim
 * @param num
 */
void printLNORM(double** matrix, int dim, int num);
/**
 *
 * @param matrix
 * @param dim
 * @param num
 */
void printJACOBI(double** matrix, int dim, int num);
/**
 *
 * @param matrix
 * @param dim
 * @param num
 * @return
 */
double** getWeightedMatrix(double** matrix, int dim, int num);
/**
 *
 * @param matrix
 * @param index1
 * @param index2
 * @param dim
 * @return
 * this function calculate the value of exp(−(||xi − xj||/2)), and returns it.
 */
double calculateWeight(double** matrix, int index1, int index2, int dim);
/**
 *
 * @param matrix
 * @param num
 * @return
 * this function calculate the value of sum(wiz) and update the Diagonal Degree Matrix d,
 * if the index i=j then dij=sumOf(wiz) for all z=1..n, otherwise dij=0.
 */
double** getDiagonalDegreeMatrix(double** matrix,int dim, int num);
/**
 *
 * @param matrix
 * @param dim
 * @return
 * Form The Normalized Graph Laplacian matrix Lnorm ∈ R^(n×n).
 * Lnorm = I − (D^(-1/2) * W * D^(-1/2))
 *       = ( (Diagonal Degree Matrix * Weighted Matrix) * Diagonal Degree Matrix )
 */
double** getLaplacianMatrix(double** matrix, int dim, int num);
/**
 *
 * @param A
 * @param n
 * @param eigenvalues
 * @return
 */
double** jacobiAlgorithm(double** A, int n, double* eigenvalues);
/**
 *
 * @param dim
 * @return
 */
double** getUnitMatrix(int dim);
/**
 *
 * @param A
 * @param pMatrix
 * @param ptMatrix
 * @param dim
 */
void RotationMat(double** A,double** pMatrix,double** ptMatrix,int dim);
/**
 *
 * @param number
 * @return
 * sign(x) = 1 if x>=0, else 0
 */
int getSign(double number);
/**
 *
 * @param matrix1
 * @param matrix2
 * @param dim
 * @return
 */
double** multiplyMatrices(double** matrix1, double** matrix2, int dim);
/**
 * @param matrix
 * @param dim
 * @return
 * Let off(A)^2 be the sum of squares of all off-diagonal elements of A respectively.
 */
double offMatrix(double** matrix, int dim);
/**
 *
 * @param matrix1
 * @param matrix2
 * @param dim
 * @return
 */
bool convergence(double** matrix1,double** matrix2,int dim);
/**
 *
 * @param V
 * @param eigenvalues
 * @param dim
 * @return
 */
double** concatenation(double **V, const double *eigenvalues, int dim);
/**
 * @param matrix
 * @param len
 * free memory function
 */
void freeMemory(double** matrix ,int len);
/**
 * @param col
 * @param row
 * @return
 * returns a new matrix with dimension col x row
 *
 */
double** createMat(int col, int row);
/**
 * @param fileName
 * @return
 * calculate the number of col in fileName.
 */
int calculateCol(char* fileName);
/**
 * @param fileName
 * @return
 * calculate the number of rows in fileName.
 */
int calculateRows(char* fileName);
/**
 * @param fileName
 * @param inputMat
 * fill 2-dimensional arrays.
 */
void fillMat(char* fileName,double** inputMat);
/**
 *
 * @param matrix
 * @param dim
 * @param num
 */
void printMatJacobi(double** matrix, int dim, int num);
/**
 *
 * @param row
 * @param col
 * @param mat
 */
void resetMat(int row,int col,double** mat);




/**
 *
 * @param array
 * @param len
 */
int getEigengapHeuristic(double* array,int len);
double** getTMatrix(double** matrix,int dim, int num, int k);
void mergeSort(double *arr_double, int len,int* arr_int);
double* slice_double(double *arr, int start, int end);
int* slice_int(int *arr, int start, int end);
void merge(double *result, double *left_double, double *right_double, int leftLen, int rightLen,int* ind,int* left_int,int* right_int);


/* kmeans algorithm in c ,printing in python*/
double diff_norm_pow2(double* vector1, double* vector2, int vector_dim);
int calc_norm(double** points, double* point, int k, int vector_dim);
double* divide(double* vector,int num,int vector_dim);
void controlPanel(int k , int max_iter, int d, int numPoints, double **all_points, double **init_centroids);

#endif