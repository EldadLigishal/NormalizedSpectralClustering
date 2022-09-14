from enum import Enum, auto
from operator import le
from xml.etree.ElementTree import tostring
import numpy as np
import pandas as pd
import sys
import myspkmeans as km


# The MAX ITER variable should be set to 300
MAXITR = 300


class Goal(Enum):
    SPK = "spk"
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"

    # check if a value exists in an enum:
    def has_value(value):
        return (value == "spk") or (value == "wam") or \
               (value == "ddg") or (value == "lnorm") or (value == "jacobi")


# k := the number of clusters required.
def execute(k, goal, input_filename):
    # reading input file
    try:
        input_matrix = pd.read_csv(input_filename, header=None)
    except:
        print("An Error Has Occurred")
        return 0

    # check if goal is valid
    if not Goal.has_value(goal):
        print("Invalid Input!")
        return 0

    matrix = input_matrix.to_numpy()
    arraylist = matrix.tolist()
    # n := number of line/rows of the input file = number of vectors = len(inputMat).
    n = len(arraylist)
    # d := number of column of an input file.
    d = len(arraylist[0])

    # Check if the data is correct
    # k must be < the number of datapoints in the file given
    if k < 0:
        print("Invalid Input!")
        return 0
    
    if Goal(goal) is Goal.SPK:
        
        if k > n:
            print("Invalid Input!")
            return 0

        Tmatrix = km.fitspkgetT(arraylist, k, d, n)
        Tmatrix_array = np.array(Tmatrix)
        k = len(Tmatrix_array[0])
        centroids = buildCentroids(k, n, Tmatrix_array)
        centroidsList = centroids.tolist()
        TmatrixList = Tmatrix_array.tolist()
        # d = k
        spk_matrix = km.fitspk(k, MAXITR, k, TmatrixList, centroidsList)
        printMatrix(spk_matrix)
    else:
        # goal == Goal.DDG or Goal.WAM or Goal.JACOBI or Goal.LNORM         
        result = km.fitope(arraylist, d, n, goal)
        printMatrix(result)


def buildCentroids(k, n, input_matrix):
    centroids = np.zeros((k, len(input_matrix[0])))
    centroids_index = np.zeros(k)
    # Select Âµ1 randomly from x1, x2, . . . , xN
    np.random.seed(0)
    random_index = np.random.choice(n, 1)
    centroids_index[0] = random_index
    centroids[0] = input_matrix[random_index]

    i = 1
    while i < k:
        D = np.zeros(n)
        # Dl = min (xl âˆ’ Âµj)^2 âˆ€j 1 â‰¤ j â‰¤ i
        for l in range(n):
            minimum = step1(i, input_matrix[l], centroids)
            D[l] = minimum

        # randomly select Âµi = xl, where P(Âµi = xl) = P(xl)
        prob = np.zeros(n)
        sum_matrix = np.sum(D)
        for j in range(n):
            prob[j] = D[j] / sum_matrix

        rand_i = np.random.choice(n, p=prob)
        centroids_index[i] = rand_i
        centroids[i] = input_matrix[rand_i]
        i += 1

    # The first line will be the indices of the observations chosen by the K-means++ algorithm
    # as the initial centroids. Observationâ€™s index is given by the first column in each input file.
    printIndex(centroids_index)
    return centroids

# vector := input_matrix[l]
# matrix := centroids
def step1(i, vector, matrix):
    min_vector = distance(vector, matrix[0])
    for j in range(0, i):
        curr_vector = distance(vector, matrix[j])
        if min_vector > curr_vector:
            min_vector = curr_vector
    return min_vector


# calculate squared distance between 2 points
def distance(vector1, vector2):
    dis = 0.0
    for i in range(len(vector1)):
        dis += ((vector1[i] - vector2[i]) * (vector1[i] - vector2[i]))
    return dis


# print 2 dimension array
def printMatrix(matrix):
    for i in range(len(matrix)):
        print(','.join([format(matrix[i][j], ".4f") for j in range(len(matrix[i]))]))


#**************************
# print 1 dimension array *
#**************************
def printIndex(array):
    for i in range(0, len(array.astype(int))):
        print(array.astype(int)[i], end="")
        if i + 1 != len(array.astype(int)):
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
    # k MUST be passed for all goals.
    execute(int(input_argv[1]), input_argv[2], input_argv[3])
else:
    print("Invalid Input!")