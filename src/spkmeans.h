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
 * @param wMatrix
 * @param dim
 * @return
 * Form the diagonal degree matrix D ∈ R^(n×n).
 * the diagonal equals to the sum of the i-th row of W.
 */
double** getDiagonalDegreeMatrix(double** wMatrix, int dim);

/**
 *
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
#endif
