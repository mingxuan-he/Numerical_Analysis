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

