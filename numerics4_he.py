import numpy as np
import numerics2_he as num2


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


def householder(vec):
    # Input: a vector x
    # Output: the householder reflector I-2P mapping x to ||x||*(1,0,0,...0) of the same length

    m = vec.size
    w = np.zeros(m)
    w[0] = num2.norm(vec,2)
    # print("input",vec)
    # print("w", w)
    v = w - vec
    # print("v", v)
    vt = np.transpose(v)
    P = np.outer(v,vt) / np.dot(vt,v)
    # print("P",P)
    I = np.identity(m)

    return I - 2 * P


def qr(matrix):
    # Input: a matrix (mxn)
    # Output: an orthogonal matrix Q and an upper triangular matrix R so that QR=matrix
    # This implementation computes the QR factorization using Householder reflectors.

    # Check Input
    matrix = np.asarray(matrix)
    m,n = np.shape(matrix)
    if m < n:
        raise ValueError("The matrix is underdetermined.")

    # Initialization
    Q = np.identity(m)
    R = matrix

    for col in range(n):
        H_hat = householder(R[col:m,col])
        H = np.block([[np.identity(col), np.zeros((col,m-col))],
                     [np.zeros((m-col,col)),H_hat]])
        R = np.matmul(H,R)
        Q = np.matmul(Q,H)

        """
        print("H_hat",H_hat)
        print("H",H)
        print("Q", Q)
        print("R", R)
        print("QR", np.dot(Q,R))
        """

    return Q, R


def qrsolve(A,b):
    # Input: (mxn) matrix A and n vector b
    # Output: x so that Ax=b if A is square or x is least squares solution if m>n
    #                resid as the error of the least squares

    # Check Inputs
    A = np.asarray(A)
    b = np.asarray(b)
    m,n = A.shape

    Q, R = qr(A)

    # For non-square matrices, extract Q_hat and R_hat
    if m != n:
        Q = Q[:,:n]
        R = R[:n,:]

    QTb = np.dot(np.transpose(Q),b)

    x = num2.backsub(R,QTb)

    resid = num2.norm(b-np.dot(A,x),2)


    return x, resid


"""
Simple Test run:

b = np.array([1,2])
bt = np.transpose(b)
print("b",b)
print("bt", bt)
print("bbt", np.outer(b,bt))
ref = householder(b)
print("reflector", ref)
print(np.dot(ref, b))

A = np.array([[1.,1],[1,2],[1,3],[1,4],[1,5]])
B = np.array([[6.,3,2,-1],[2,1,-2,1],[2,9,4,2],[7,3,-1,1],[6,-6,8,3]])
print(A)
qr(A)
"""