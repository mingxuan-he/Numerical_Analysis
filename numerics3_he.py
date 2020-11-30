import numpy as np
import numerics0_he as num0
import numerics4_he as num4


def newtondd(xdata,ydata):
    # xdata, ydata: vectors that need to be interploated as (x[i],y[i]).

    # Check inputs
    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)
    x_obs = len(xdata)
    y_obs = len(ydata)
    if x_obs != y_obs:
        raise ValueError("Data points do not match.")

    # Initialization
    ddtable = np.zeros((x_obs,x_obs))
    coeffs = np.zeros(x_obs)
    ddtable[:, 0] = ydata
    coeffs[0] = ddtable[0,0]

    # Compute divided differences
    for col in range(1,x_obs):
        for row in range(col,x_obs):
            ddtable[row,col] = (ddtable[row,col-1] - ddtable[row-1,col-1]) / (xdata[row] - xdata[row - col])
        coeffs[col] = ddtable[col,col]

    # print(ddtable)
    return coeffs


def newtonInterp(xdata,ydata):
    # xdata, ydata: vectors that need to be interploated as (x[i],y[i]).
    # Output: poly which is a function to evaluate the minimal degree interpolating polynomial defined by this data set.

    coeff = newtondd(xdata,ydata)
    poly = lambda x: num0.polynest(x,coeff,xdata)
    return poly


def chebyshevRoots(num_roots):
    # num_roots: number of Chebyshev roots needed

    roots = np.zeros(num_roots)

    for i in range(num_roots):
        roots[i] = np.cos((2*i+1)*np.pi / (2*num_roots))

    return roots


def chebyshevInterp(function, num_roots, interval=[-1,1], interp_method="Newton"):
    # function: a function of one variable
    # num_roots: number of Chebyshev roots needed
    # interval: the interpolating interval for x, default is [-1,1]

    # Calculate nodes on [-1,1]
    roots = chebyshevRoots(num_roots)
    # print(roots)

    # Change roots according to interval if necessary
    if interval != [-1,1]:
        [a,b] = interval
        roots = (b-a)/2 * roots + (b+a)/2
        # print(roots)

    # Take sample values from the original function
    yval = function(roots)

    # Interpolate using the given method
    # For now it's just Newton, but We can add more interpolate methods later
    if interp_method == "Newton":
        chebyshevfunc = newtonInterp(roots, yval)
    else:
        raise ValueError("Choose an appropriate interpolation method from the following: Newton.")

    return chebyshevfunc

"""
Simple test run:
def func(x):
    return np.sin(x)

chebyshevInterp(func,6,interval=[0,2*np.pi])
"""


def cubiccoeff(xdata, ydata, end_condition="natural", df = lambda x: 0):
    # xdata, ydata: vectors that need to be interploated as (x[i],y[i]). needs to be sorted
    # end_condition is the specified endpoint condition, default is natural
    # df is the derivative of the interpolation function, default is f': R->{0}

    # Check inputs
    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)
    m = len(xdata)
    n = len(ydata)
    if m !=n:
        raise ValueError("Data points do not match.")

    # Initialization
    left_endpoint = xdata[0]
    right_endpoint = xdata[-1]
    coeff = np.zeros((n-1,4))
    x_dist = np.zeros(n-1)
    y_dist = np.zeros(n-1)
    a = b = d = np.zeros(n-1)

    # Calculate length of intervals (delta)
    for i in range(0,n-1):
        x_dist[i] = xdata[i+1] - xdata[i]
        y_dist[i] = ydata[i+1] - ydata[i]

    # Calculate coefficients a,b,c,d

    # a
    a = ydata[:-1]

    # c
    if end_condition == "natural":

        # Make the delta matrix
        delta_mat = np.identity(n)
        for row in range(0,n-2):
            delta_mat[row+1,row] = x_dist[row]
            delta_mat[row+1,row+1] = 2 * (x_dist[row] + x_dist[row+1])
            delta_mat[row+1,row+2] = x_dist[row+1]
        #print("matrix",delta_mat)

        # Make the delta vector (as product)
        delta_vec = np.zeros(n)
        delta_vec[0] = 0
        delta_vec[-1] = 0
        for i in range(0,n-2):
            delta_vec[i+1] = 3 * (y_dist[i+1]/x_dist[i+1] - y_dist[i]/x_dist[i])
        #print("vec",delta_vec)
        c,resid = num4.qrsolve(delta_mat,delta_vec)

    elif end_condition == "clamped":

        if df == (lambda x: 0):
            raise ValueError("Derivative is required to use clamped cubic spline.")

        # Make the delta matrix
        delta_mat = np.identity(n)
        delta_mat[0,0] = 2 * x_dist[0]
        delta_mat[0,1] = x_dist[0]
        for row in range(0,n-2):
            delta_mat[row + 1, row] = x_dist[row]
            delta_mat[row + 1, row + 1] = 2 * (x_dist[row] + x_dist[row + 1])
            delta_mat[row + 1, row + 2] = x_dist[row + 1]
        delta_mat[-1,-2] = x_dist[-1]
        delta_mat[-1,-1] = 2 * x_dist[-1]

        # Make the delta vector (as product)
        delta_vec = np.zeros(n)
        delta_vec[0] = 3 * (y_dist[0] / x_dist[0] - df(left_endpoint))
        delta_vec[-1] = 3 * (df(right_endpoint) - y_dist[-1] / x_dist[-1])
        for i in range(0, n - 2):
            delta_vec[i + 1] = 3 * (y_dist[i + 1] / x_dist[i + 1] - y_dist[i] / x_dist[i])
        # print("vec",delta_vec)
        c, resid = num4.qrsolve(delta_mat, delta_vec)

    else:
        raise ValueError("Choose an appropriate endpoint condition from the following: natural, clamped.")

    # d
    for i in range(n-1):
        d[i] = (c[i+1] - c[i]) / (3 * x_dist[i])

    # b
    for i in range(n-1):
        b[i] = y_dist[i] / x_dist[i] - (x_dist[i] / 3) * (2 * c[i] + c[i+1])

    coeff = np.stack((a,b,c[:-1],d))
    """
    print("a",a)
    print("b",b)
    print("c",c)
    print("d",d)
    """

    return coeff


def cubic_eval(xdata, ydata, x, end_condition="natural", df = lambda x: 0):
    # xdata, ydata, end_condition, df: same as cubiccoeff
    # x: a point for evaluation

    # Binary search for the position of x
    left_idx = 0
    right_idx = len(xdata) - 1
    mid_idx = 0

    while xdata[left_idx+1] < x or xdata[right_idx-1] > x:
        mid_idx = (right_idx - left_idx) // 2

        if xdata[mid_idx] < x:
            left_idx = mid_idx

        elif xdata[mid_idx] > x:
            right_idx = mid_idx

        # x is a given point
        else:
            return ydata[mid_idx]
        #print(left_idx, right_idx)

    # Evaluate using nested polynomial
    coeff = cubiccoeff(xdata,ydata,end_condition,df)
    polycoeff = coeff[:,left_idx]
    evaluation = num0.polynest(x,polycoeff,xdata)

    return evaluation
