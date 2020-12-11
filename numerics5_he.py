import numpy as np


def simpson(function, a, b):
    # function: a function to integrate on the interval [a,b]
    # a,b: left and right endpoint of the integrating interval

    # Check inputs
    a = float(a)
    b = float(b)
    if a >= b:
        raise ValueError("[a,b] must be a real interval.")

    # Compute and evaluate at endpoints and midpoint
    h = (b-a) / 2
    x = [a, a + h, b]
    y = function(x)

    # Compute definite integral
    integral = (h / 3) * (y[0] + 4 * y[1] + y[2])

    return integral


def compositeSimpson(function, n, a, b):
    # function: a function to integrate on the interval [a,b]
    # n: number of intervals to partition the interval
    # a,b: left and right endpoint of the integrating interval

    # Check inputs
    n = int(n) - (int(n) % 2)
    a = float(a)
    b = float(b)
    if a >= b:
        raise ValueError("[a,b] must be a real interval.")

    # Compute and evaluate at sample points
    x = np.linspace(a,b,n+1)
    #print("Sample points", x)
    h = x[1] - x[0]
    y = function(x)

    # Compute definite integral
    integral = 0
    # endpoints
    integral += y[0] + y[-1]
    # x_{2i}
    for i in range(2,n-1,2):
        integral += 2*y[i]
    # x_{2i+1}
    for i in range(1,n,2):
        integral += 4*y[i]
    integral = integral * (h/3)

    return integral


def adaptiveSimpson(function, a, b, tolerance=1.0e-12, Sab=0.0, recursion_counter=0, max_depth=20):
    # function: a function to integrate on the interval [a,b]
    # a,b: left and right endpoint of the integrating interval
    # tolerance: tolerance distributed to the specific sub-interval
    # Sab: Simpson integration value on a current interval [a_i,b_i]
    # recursion_counter: counts the current recursion level

    # Check inputs
    if recursion_counter == 0:
        a = float(a)
        b = float(b)
        if a >= b:
            raise ValueError("[a,b] must be a real interval.")
    elif recursion_counter > max_depth:
        return Sab

    # Midpoint
    c = (a+b)/2

    # Compute simpson's on left and right interval, then sum
    left_approx = simpson(function,a,c)
    right_approx = simpson(function,c,b)
    sum_approx = left_approx + right_approx

    """
    print("depth: %d" % recursion_counter)
    print("Sab", Sab)
    print("sum_approx", sum_approx)
    print("tolerance", tolerance)
    """
    # Compare the left-right sum with previous approximation
    # If error is small, done.
    if abs(Sab - sum_approx) <= 10 * tolerance:
        #print("we're good")
        return sum_approx
    else:
        # If error is large, do next-level recursion and add them up:
        # interval become [a,c] and [c,b], Sab become the computed approximation on [a,c] and [c,b];
        # tolerance become half the current tolerance, recursion counter ++
        #print("going next")
        return adaptiveSimpson(function,a,c,tolerance=tolerance/2,Sab=left_approx,recursion_counter=recursion_counter+1) + adaptiveSimpson(function,c,b,tolerance=tolerance/2,Sab=right_approx,recursion_counter=recursion_counter+1)


# for legendre polynomials: np.polynomial.legendre.leggauss(n)

def gaussQuad(function,n,a,b):
    # function: a function to integrate on the interval [a,b]
    # a,b: left and right endpoint of the integrating interval
    # n: number of sample points taken

    # Check inputs
    a = float(a)
    b = float(b)
    if a >= b:
        raise ValueError("[a,b] must be a real interval.")

    # change interval
    if [a,b] != [-1,1]:
        function = lambda t: function(t * (b-a)/2 + (b+a)/2)

    nodes,coeffs = np.polynomial.legendre.leggauss(n)

    integral = 0
    for i in range(n):
        integral += coeffs[i] * function(nodes[i])

    return integral

