import numpy as np
import matplotlib.pyplot as plt


def bisection(function, left_endpoint, right_endpoint, tolerance=1.0e-6, maxIter=100):
    # function is a numerical function with single variable/input x
    # left endpoint (a) and right endpoint (b) are both numbers on the real line
    # tolerance is the error/precision we are looking for, default is 1.0e-6
    # maxIter is the max number of iterations, stops the algorithm once reached
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
    # function is a numerical function with single variable/input x
    # x init is the initial input value for x (guessed)
    # tolerance is the error/precision we are looking for, default is 1.0e-6
    # maxIter is the max number of iterations, stops the algorithm once reached
    # Important note: the method might not find a convergent fixed point if the conditions of Thm 1.2.a are not met.
    # If not, first convert the function into g(x) which FPI can be applied

    # Initialization
    x = xinit
    fx = function(x)
    iter_count = 0
    err = abs(x)
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
Simple test run:

def f(x):
    return pow(x,2)

print(fixedpt(f,0.1))
"""


