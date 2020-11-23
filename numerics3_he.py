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

    print(ddtable)
    return coeffs


def newtonInterp(xdata,ydata):
    # Input: xdata and ydata vectors that need to be interploated as (x[i],y[i]).
    # Output: poly which is a function to evaluate the minimal degree interpolating polynomial defined by this data set.

    coeff = newtondd(xdata,ydata)
    poly = lambda x: num0.polynest(x,coeff,xdata)
    return poly

