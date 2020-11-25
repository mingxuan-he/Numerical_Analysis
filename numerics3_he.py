import numpy as np
import numerics0_he as num0


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
