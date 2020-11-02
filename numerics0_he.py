import numpy as np

def polynest (x, a, b=[]):
    # x is the input to the polynomial
    # a is the vector of n+1 coefficients of the polynomial
    # b is the vector of n base points. Default is empty
    # This function evaluate a polynomial with coefficients a and basepoints b at the point x and returns that value as y

    # Read a,b as numpy arrays
    if b is None:
        b = []
    a = np.asarray(a)
    b = np.asarray(b)

    # Number of dimensions on a,b
    adims = a.ndim
    bdims = b.ndim

    # Input error check
    if adims != 1:
        raise ValueError("The coefficients (2nd input) must be passed as a vector (1 dimension).")
    if bdims != 1:
        raise ValueError("The base points (optional 3nd input) must be passed as a vector (1 dimension).")

    n = a.size
    numbasepoints = b.size
    y = a[n-1]

    # Execute algorithm
    if numbasepoints:   # If there are base points
        for idx in range(n-2, -1, -1):
            y = y * (x-b[idx])+a[idx]

    else:   # If there are no base points:
        for idx in range(n-2,-1,-1):
            y = y * x + a[idx]

    return y



