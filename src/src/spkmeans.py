import enum
import math
import numpy as np
import pandas as pd
import sys
import mykmeanssp as km
EPSILON = 10**-10   # good epsilon?


class Goal(enum):
    SPK = "spk"
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"

    def has_value(self, value):
        return (value == self.SPK) or (value == self.WAM) or (
            value == self.DDG) or (self.LNORM) or (self.JACOBI)



# k := the number of clusters required.
def execute(k, goal,input_filename):

    input_matrix = pd.read_csv(input_filename)
    
    # check if goal is valid
    if not Goal.has_value(goal):
        print("Invalid Input!")
        return 0

    matrix = input_matrix.to_numpy()
    # n := number of rows of the input file.
    # d := number of columns of the input file.
    n = matrix.shape(0)
    d = matrix.shape(1)

    # Check if the data is correct
    if (k >= n) or (k < 0):
        print("Invalid Input!")
        return 0


    if k == 0:
        k = # getEigengapHeuristic
    
    
# last editing
# ---------------------------------------------------------------------------------------------

    # centroids Âµ1, Âµ2, ... , ÂµK âˆˆ R^d where 1<K<N.
    centroids = buildCentroids(k, n, input_matrix)

    matrix = km.fit(k, 300, EPSILON, n, d, input_array, centroids.tolist())
    printMatrix(np.array(matrix))


def buildCentroids(k, n, input_matrix):
    centroids = np.zeros((k, len(input_matrix[0])))
    centroids_index = np.zeros(k)
    # Select Âµ1 randomly from x1, x2, . . . , xN
    np.random.seed(0)
    random_index = np.random.choice(n, 1)
    # ERROR: setting an array element with a sequence.
    centroids_index[0] = random_index
    centroids[0] = input_matrix[random_index]

    i = 1
    while i < k:
        d = np.zeros(n)
        # Dl = min (xl âˆ’ Âµj)^2 âˆ€j 1 â‰¤ j â‰¤ i
        for _ in range(n):
            d[_] = step1(i, input_matrix[_], centroids)
        # randomly select Âµi = xl, where P(Âµi = xl) = P(xl)
        prob = np.zeros(n)
        sum_matrix = np.sum(d)
        for j in range(n):
            prob[j] = d[j]/sum_matrix
       #     prob[j] = step2(d[j], sum_matrix)
        rand_i = np.random.choice(n, p=prob)
        centroids_index[i] = rand_i
        centroids[i] = input_matrix[rand_i]
        i += 1

    # The first line will be the indices of the observations chosen by the K-means++ algorithm
    # as the initial centroids. Observationâ€™s index is given by the first column in each input file.
    printIndex(centroids_index)
    return centroids
 

def step1(i, vector, matrix):
    min_vector = math.pow(np.linalg.norm(np.subtract(vector, matrix[0])), 2)
    for j in range(i):
        curr_vector = np.linalg.norm(np.subtract(vector, matrix[j]))
        curr_vector = np.power(curr_vector, 2)
        if min_vector > curr_vector:
            min_vector = curr_vector
    return min_vector


def step2(vector, sum_matrix):
    result = np.divide(vector, sum_matrix)
    return result


def printMatrix(arr):
    for i in range(0, len(arr)):
        for j in range(0, len(arr[0])):
            print(np.round(arr[i][j], 4), end="")
            if j + 1 != len(arr[0]):
                print(",", end="")
        print()


def printIndex(matrix):
    for i in range(0, len(matrix.astype(int))):
        print(matrix.astype(int)[i], end="")
        if i + 1 != len(matrix.astype(int)):
            print(",", end="")
    print()


# main
# sys.argv is the list of command-line arguments.
# len(sys.argv) is the number of command-line arguments.
# sys.argv[0] is the program i.e. the script name.
input_argv = sys.argv
input_argc = len(sys.argv)
if input_argc == 4:
    # (k, goal, inputFile1)
    execute(int(input_argv[1]), input_argv[2], input_argv[3])
else:
    # (k, epsilon, inputFile1, inputFile2)
    print("Invalid Input!")