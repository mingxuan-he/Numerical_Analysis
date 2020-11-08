import numpy as np


# Row operations:
def rowdiff(matrix,k,j,scale):
    # matrix: a 2 dimensional numpy array
    # k,j: row numbers
    # scale: a constant to scale row j
    # replace row k by (k - scale * j)
    matrix[k,:] = matrix[k,:] - scale * matrix[j,:]
    return matrix


def rowswap(matrix,k,j):
    # matrix: a 2 dimensional numpy array
    # k,j: row numbers
    matrix[[k,j],:] = matrix[[j,k],:]
    return matrix


def rowscaling(matrix,k,scale):
    # matrix: a 2 dimensional numpy array
    # k: row number
    # scale: a constant to scale row k
    matrix[k,:] = scale * matrix[k,:]
    return matrix




