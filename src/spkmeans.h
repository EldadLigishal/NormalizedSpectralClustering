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
 * @param n
 * @return
 * this func Form the weighted adjacency matrix W ∈ R^(n×n),
 * for the matrix given as input and returns it.
 */
double** getWeightedMatrix(double** matrix, int n);

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

#endif
