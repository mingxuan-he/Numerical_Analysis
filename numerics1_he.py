import numpy as np


def bisection(function, left_endpoint, right_endpoint, tolerance=1.0e-6, maxIter=100):
    # function: a numerical function with single variable/input x
    # left endpoint (a) and right endpoint (b): both numbers on the real line
    # tolerance: the error/precision we are looking for, default is 1.0e-6
    # maxIter: the max number of iterations, stops the algorithm once reached
    # (a,b) is a real interval, meaning the inputs must satisfy that a<b

    # Check if input is a proper interval
    if left_endpoint >= right_endpoint:
        raise ValueError("The left endpoint must be to the left of the right endpoint.")

    # Evaluate the function at the initial endpoints, get sign
    sign_fa = np.sign(function(left_endpoint))
    sign_fb = np.sign(function(right_endpoint))

    # Check if the input interval contains a root
    if sign_fa == 0: return left_endpoint
    if sign_fb == 0: return right_endpoint
    if sign_fa * sign_fb == 1:
        raise ValueError("The given interval might not contain a root.")

    # Initialization
    iter_count = 0
    midpoint = 0
    midpoints = np.zeros((maxIter, 1))

#   while function(left_endpoint) * function(right_endpoint) <0: # Under this condition, the function would loop infinitely
    while (right_endpoint - left_endpoint) > tolerance and iter_count <= maxIter:
        midpoint = (left_endpoint + right_endpoint) * 0.5
        # Record the new midpoint as an array entry
        midpoints[iter_count] = midpoint
        # If the midpoint is a root, return it
        if function(midpoint) == 0:
            return midpoint
        # If left half contains a root, select the left half
        elif function(left_endpoint) * function(midpoint) < 0:
            right_endpoint = midpoint
        # Otherwise, select the right half
        else:
            left_endpoint = midpoint
        iter_count += 1

        # print(left_endpoint, right_endpoint)

    # Return both the final root and the iterated midpoints
    return midpoint, midpoints[0:iter_count]


"""
Bisection simple test run:

def f(x):
    return pow(x,2) - 1

print(bisection(f,0,2))
"""


def fixedpt(function, xinit, tolerance=1.0e-6, maxIter=100):
    # function: a numerical function with single variable/input x
    # xinit: the initial input value for x (guessed)
    # tolerance: the error/precision we are looking for, default is 1.0e-6
    # maxIter: the max number of iterations, stops the algorithm once reached
    # Important note: The conditions in Thm 1.2.a must be satisfied.
    # If the original function does not satisfy the conditions, first convert the function into g(x) which FPI can be applied.

    # Initialization
    x = xinit
    fx = function(x)
    iter_count = 0
    err = abs(fx)
    roots = np.zeros((maxIter,1))
    roots[0] = xinit

    while err > tolerance and iter_count < maxIter:
        # update x, f(x) iteratively, then compute the error
        x = fx
        fx = function(fx)
        err = abs(fx - x)
        roots[iter_count+1] = x
        iter_count += 1
        #print(x)
        #print(fx)

    return x, roots[0:iter_count+1]




"""
Newton simple test run:

def f(x):
    return pow(x,2)

print(fixedpt(f,0.1))
"""


def newton(function, dfunction, xinit, tolerance=1.0e-6, maxIter=100):
    # function: a numerical function with single variable/input x
    # xinit: the initial input value for x (guessed)
    # dfunction: the derivative of function
    # tolerance: the error/precision we are looking for, default is 1.0e-6
    # maxIter: the max number of iterations, stops the algorithm once reached
    # Important note: the derivative at the target root cannot be zero (simple root) so that the function converges locally and quadratically

    # Initialization
    x = xinit
    fx = function(x)
    dfx = dfunction(x)
    iter_count = 0
    err = fx / dfx
    roots = np.zeros((maxIter, 1))
    roots[0] = xinit

    while err > tolerance and iter_count < maxIter:
        # update x, f(x) iteratively, then compute the error
        x = x - (fx / dfx)
        fx = function(x)
        dfx = dfunction(x)
        err = abs(fx / dfx)
        roots[iter_count+1] = x
        iter_count += 1
        #print(x)
        #print(fx)

    return x, roots[0:iter_count+1]

"""
Secant simple test run:

def f(x):
    return x**3 -1

def df(x):
    return 3 * x**2

print(newton(f,df,2))
"""


def secant(function, xinit, tolerance=1.0e-6, maxIter=100):
    # function: a numerical function with single variable/input x
    # xinit: the initial input value for x (guessed)
    # tolerance: the error/precision we are looking for, default is 1.0e-6
    # maxIter: the max number of iterations, stops the algorithm once reached
    # Important note: the secant method is similar to the Newton's method except that it does not require derivative of f

    # Initialization
    x = xinit
    fx = function(x)
    sec = fx
    iter_count = 0
    err = fx / sec
    roots = np.zeros((maxIter, 1))
    roots[0] = xinit
    fvals = np.zeros((maxIter, 1))
    fvals[0]= fx

    while err > tolerance and iter_count < maxIter:
        # update x, f(x) iteratively, then compute the error
        x = x - (fx / sec)
        fx = function(x)
        sec = (fx - fvals[iter_count]) / (x - roots[iter_count])
        err = abs(fx / sec)
        roots[iter_count+1] = x
        fvals[iter_count+1] = fx
        iter_count += 1
        print(x)
        print(fx)

    return x, roots[0:iter_count+1]


"""
Simple test run:

def f(x):
    return x**3 -1
    
print(secant(f,4))
"""