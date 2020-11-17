import numpy as np


def leastSquares_lu(A,b):

    A = np.asarray(A)
    #m,n = A.shape

    At = np.transpose(A)
    AtA = np.dot(At, A)
    y = np.dot(At, b)
    print("A^tA =")
    print(AtA)
    print("A^tb =")
    print(y)

    x = num2.lusolve(AtA,y)

    r = b - np.dot(A,x)

    return x,r


