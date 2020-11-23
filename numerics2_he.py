import numpy as np


def norm(vec, l):
    # vec: a vector (1D numpy array)
    # l: dimension of the norm, l = -1 is the l-infinity norm
    # Returns the l-norm of a vector (a scalar)

    # Check input
    vec = np.asarray(vec)

    if l == -1:
        return np.amax(abs(vec))

    elif l > 0:
        powerf = lambda t: t**l
        return np.sum(powerf(abs(vec))) ** (1 / l)


# Row operations:
def rowdiff(matrix,k,j,scale=1.0):
    # matrix: a 2 dimensional numpy array
    # k,j: row numbers
    # scale: a constant to scale row j
    # replace row k by (k - scale * j)
    matrix[k,:] = matrix[k,:] - scale * matrix[j,:]
    return matrix


def rowswap(matrix,k,j):
    # matrix: a 2 dimensional numpy array
    # k,j: row numbers
    # swap 2 rows in a matrix
    matrix[[k,j],:] = matrix[[j,k],:]
    return matrix


def rowscale(matrix,k,scale=1.0):
    # matrix: a 2 dimensional numpy array
    # k: row number
    # scale: a constant to scale row k
    # scale a row in a matrix
    matrix[k,:] = scale * matrix[k,:]
    return matrix


def noswapLU(matrix):
    # matrix: an n by n matrix
    # This function performs LU decomposition on the matrix
    # This function does not use and swap operation

    # Check input matrix
    A = np.asarray(matrix)
    (m,n) = A.shape
    if (m != n):
        raise ValueError("The input matrix must be square. Your matrix has size %d x %d." % (m,n))

    # initialize square matrix L
    L = np.zeros((m,n))

    for col in range(n):
        L[col,col] = 1.0    # first write 1 in the diagonal entry of col in L
        pivot = np.float(A[col,col]) # compute the pivot of the row combination
        # for each row index greater than the colum index, do row combination with this scaling factor
        for row in range(col+1,n):
            scalefactor = A[row,col] / pivot
            L[row,col] = scalefactor
            A = rowdiff(A,row,col,scalefactor)

    return L, A


def backsub(U,b):
    # This function takes in an upper triangular matrix U and solves the system Ux=b using backward substitution
    # Returns the solution vector x

    # Check inputs
    U = np.asarray(U)
    b = np.asarray(b)
    n = b.size
    (mU,nU) = U.shape
    if (mU != nU) or (n != nU):
        raise ValueError("The dimensions of the input are incorrect. U must be square and b must be a vector. The number of columns of U must match the length of b.")

    # Initialize the output solution vector
    x = np.zeros(b.shape)

    # Read the last entry
    x[n-1] = b[n-1] / np.float(U[n-1,n-1])
    # Then work backwards from row n-2 up to row 0
    for row in range(n-2,-1,-1):
        # We have found the values of x with index greater than the current row
        # The upper triangular matrix may have factors for those values of x
        # we complete the short inner product, subtract it from the associated value in b, and divide by the coefficient in U
        x[row] = (b[row] - np.dot(U[row, row+1:n], x[row+1:n])) / U[row, row]

    return x


def forwardsub(L,b):
    # This function takes in an lower triangular matrix L and solves the system Lx=b using forward substitution
    # Returns the solution vector x

    # Check inputs
    L = np.asarray(L)
    b = np.asarray(b)
    n = b.size
    (mL, nL) = L.shape
    if (mL != nL) or (n != nL):
        raise ValueError("The dimensions of the input are incorrect. U must be square and b must be a vector. The number of columns of L must match the length of b.")

    # Initialize the output solution vector
    x = np.zeros(b.shape)

    # Read the first entry
    x[0] = b[0] / np.float(L[0, 0])
    # Then work forward from row 1 up to row n-1
    for row in range(1,n):
        # We have found the values of x with index smaller than the current row
        # The lower triangular matrix may have factors for those values of x
        # we complete the short inner product, subtract it from the associated value in b, and divide by the coefficient in U
        x[row] = (b[row] - np.dot(L[row, 0:row], x[0:row])) / L[row, row]

    return x


def fbsolve(L,U,b):
    # L: a lower triangular matrix
    # U: an upper triangular matrix
    # b: a vector
    # This function conducts a forward backward solve of the system LUx=b
    # Returns the solution vector x

    # Check inputs
    L = np.asarray(L)
    U = np.asarray(U)
    b = np.asarray(b)
    n = b.size
    (mU,nU) = U.shape
    (mL,nL) = L.shape
    if (mL != nL) or (n != nL) or (mU != nU) or (nU != n):
        raise ValueError("The dimensions of the input are incorrect. U and L must be square and b must be a vector. The number of columns in U and L must match the length of b.")

    y = forwardsub(L,b)
    x = backsub(U,y)

    return x


def noswaplusolve(A,b):
    # Solve the system Ax=b using LU factorization
    # Returns the solution vector x

    L,U = noswapLU(A)
    x = fbsolve(L,U,b)

    return x


def lu(matrix):
    # matrix: an n by n matrix
    # This function performs PA=LU decomposition on the matrix
    # Returns L, U, P

    # Check input:
    A = np.asarray(matrix)
    (m,n)= A.shape
    if m != n:
        raise ValueError("The input matrix must be square. Your matrix has size %d x %d" % (m,n))

    # Initialize the matrices L and P
    L = np.zeros((m,n))
    P = np.identity(n)

    # Move through n-1 columns
    for col in range(n):
        # find the row with the largest magnitude entry in this column on or below the diagonal
        pivot_index = np.argmax(np.absolute(A[col:n,col]))
        # since p is the index in that diagonal or below, we need to shift it by col-1
        pivot_index += col
        # we need the value of this absolute argmax as our pivot
        pivot = np.float(A[pivot_index,col])

        # if pivot==0, the matrix is singular; if not, we row swap to get the preferred pivot
        if pivot:
            if pivot_index != col:
                # swap A, the permutation matrix P, and the current version of L
                A = rowswap(A, pivot_index, col)
                L = rowswap(L, pivot_index, col)
                P = rowswap(P, pivot_index, col)
        else:
            raise ValueError("The matrix provided is singular so the decomposition fails.")

        for row in range(col+1,n):
            scalefactor = A[row,col] / pivot
            L[row,col] = scalefactor
            A = rowdiff(A,row,col,scalefactor)
            #print(A)

    # because of the row swapping, we need to add the ones along diagonal at the end
    L += np.identity(n)

    return L,A,P


def lusolve(A,b):
    # solves an equation Ax = b and return x by first executing an PA=LU factorization
    # and then running forward and backward substitution

    L,U,P = lu(A)
    x = fbsolve(L,U,np.dot(P,b))

    return x


def add_decomp(A):
    # Performs additive decomposition on a square matrix A
    # Returns a lower triangular matrix L, an upper triangular matrix U, and a diagonal matrix D.

    # Check input:
    A = np.asarray(A)
    (m,n)= A.shape
    if m != n:
        raise ValueError("The input matrix must be square. Your matrix has size %d x %d" % (m,n))

    # Initialize matrices L,U,D
    L = np.zeros((m, n))
    U = np.zeros((m, n))
    D = np.zeros((m, n))

    # Performs additive decomposition
    for i in range(n):
        D[i, i] = A[i, i]

    for i in range(n):
        for j in range(i):
            L[i, j] = A[i, j]

    for i in range(n):
        for j in range(i + 1, n):
            U[i, j] = A[i, j]

    """
    # Print matrices
    print(L)
    print(U)
    print(D)
    """

    return L,U,D


def jacobi(A, b, xinit, tolerance=1.0e-6, maxIter=100):
    # A: an n by n matrix
    # b: a vector
    # xinit: the initial input value for x (guessed)
    # tolerance: the error/precision we are looking for, default is 1.0e-6
    # maxIter: the max number of iterations, stops the algorithm once reached
    # This function solves Ax=b iteratively using Jacobi Method (additive decomposition), returns the solution vector x

    # Check input:
    A = np.asarray(A)
    (m,n)= A.shape
    if m != n:
        raise ValueError("The input matrix must be square. Your matrix has size %d x %d" % (m,n))
    b = np.asarray(b)
    xinit = np.asarray(xinit)

    # Initialize variables
    x = xinit
    iter_count = 0
    roots = np.zeros((maxIter, n))
    roots[0] = xinit
    err = norm((b - np.dot(A, x)), -1)

    while err > tolerance and iter_count < maxIter-1:

        for i in range(n): # each entry of x
            T_sum = 0
            for j in range(n):
                if j != i:
                    T_sum += A[i,j] * roots[iter_count][j]
                    #print("For i = %d, j = %d, T sum is %f" % (i,j,T_sum))

            x[i] = (b[i] - T_sum) / A[i,i]
        print("x", x)
        roots[iter_count+1] = x
        err = norm(b - np.dot(A, x), -1)
        #print(err)
        iter_count += 1

    return x


def gausssiedel(A, b, xinit, tolerance=1.0e-6, maxIter=100):
    # A: an n by n matrix
    # b: a vector
    # xinit: the initial input value for x (guessed)
    # tolerance: the error/precision we are looking for, default is 1.0e-6
    # maxIter: the max number of iterations, stops the algorithm once reached
    # This function solves Ax=b iteratively using Gauss-Siedel Method (additive decomposition), returns the solution vector x

    # Check input:
    A = np.asarray(A)
    (m,n) = A.shape
    if m != n:
        raise ValueError("The input matrix must be square. Your matrix has size %d x %d" % (m,n))
    b = np.asarray(b)
    xinit = np.asarray(xinit)

    # Initialize variables
    x = xinit
    iter_count = 0
    roots = np.zeros((maxIter, n))
    roots[0] = xinit
    err = norm((b - np.dot(A, x)), -1)

    while err > tolerance and iter_count < maxIter-1:

        for i in range(n): # each entry of x
            T_sum = 0
            for j in range(i):
                T_sum += A[i,j] * x[j]
            for j in range(i+1,n):
                T_sum += A[i,j] * roots[iter_count][j]

            x[i] = (b[i] - T_sum) / A[i,i]

        #print(x)
        roots[iter_count+1] = x
        err = norm((b - np.dot(A, x)), -1)
        iter_count += 1

    return x
