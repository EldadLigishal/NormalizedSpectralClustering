#ifndef spkmeans_h
#define spkmeans_h

typedef enum {SPK,WAM, DDG, LNORM, JACOBI
} Goal;

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc,char **argv);

/**
 *
 * @param matrix
 * @param dim
 * @return
 * this func Form the weighted adjacency matrix W ∈ R^(n×n),
 * for the matrix given as input and returns it.
 */
double** getWeightedMatrix(double** matrix, int dim);
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
 * @param dim
 * @return
 * this function calculate the value of sum(wiz) and update the Diagonal Degree Matrix d,
 * if the index i=j then dij=sumOf(wiz) for all z=1..n, otherwise dij=0.
 */
double** getDiagonalMatrix(double** matrix, int dim);
/**
 *
 * @param matrix
 * @param dim
 * @return
 * Form the diagonal degree matrix D ∈ R^(n×n).
 * the diagonal equals to the sum of the i-th row of W.
 */
double** getDiagonalDegreeMatrix(double** matrix, int dim);
 /**
  *
  * @param matrix
  * @param dim
  * @return
  * Form The Normalized Graph Laplacian matrix Lnorm ∈ R^(n×n).
  * Lnorm = I − (D^(-1/2) * W * D^(-1/2))
  *       = ( (Diagonal Degree Matrix * Weighted Matrix) * Diagonal Degree Matrix )
  */
double** getLaplacianMatrix(double** matrix, int dim);
/**
 *
 * @param matrix1
 * @param matrix2
 * @param dim
 * @return
 * multiplying matrices dimXdim.
 * the number of columns in the first matrix must be equal
 * to the number of rows in the second matrix.
 */
double** multiplyMatrices(double** matrix1, double** matrix2, int dim);
/**
 *
 * @param number
 * @return
 * sign(x) = 1 if x>=0, else 0
 */
int getSign(double number)
/**
 * @param matrix
 * @param dim
 * @return
 * The transpose of a matrix is found by interchanging its rows into columns.
 */
double** transpose(double** matrix, int dim);
/**
 * @param matrix
 * @param dim
 * @return
 * Let off(A)^2 be the sum of squares of all off-diagonal elements of A respectively.
 */
double offOfMat(double** matrix, int dim);
/**
 * @param array
 * @param len
 * @return
 * The Eigengap Heuristic
 * In order to determine the number of clusters k, we will use eigengap heuristic.
 */
int getEigengapHeuristic(double* array,int len);
/**
 * @param matrix
 * @param dim
 * @return
 * checks if a matrix is diagonal.
 * diagonal matrix is a square matrix that consists of all zeros off the main diagonal.
 */
bool isDiagonalMatrix(double** matrix, int dim);

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

/*****************/
void printMat(double** matrix, int dim);

#endif
