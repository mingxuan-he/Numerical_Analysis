import numpy as np
import numerics_he as num
import matplotlib.pyplot as plt


# Successive parabolic interpolation
def SPI(function, x_inits, tolerance=1.0e-6, maxIter=100, plot=False, plot_range=np.arange(0,10,0.01)):
    # function: the objective function
    # x_inits: three initial approximations
    # tolerance: maximum error allowed
    # maxIter: maximum number of iterations
    # plot: whether to plot the parabolas
    # plot_range: if plot=True, the input range for the plot

    # Check Inputs
    x_inits = np.asarray(x_inits)
    if x_inits.size != 3:
        raise ValueError("Input three initial approximations.")

    # Initialization
    k = 0
    x = np.zeros(maxIter+4)
    x[0:3] = x_inits
    y = np.zeros(maxIter+4)
    y[0:3] = function(x_inits)
    err = abs(y[2]-y[1])

    while err > tolerance and k <= maxIter:
        #print('Iteration',k)
        r = x[k]
        s = x[k+1]
        t = x[k+2]
        fr = y[k]
        fs = y[k+1]
        ft = y[k+2]

        if plot:
            # Plotting parabolas
            xdata = np.array([r,s,t])
            ydata = np.array([fr,fs,ft])
            parabola = num.newtonInterp(xdata,ydata)

            plt.title("%dth iteration" % k)
            plt.plot(xdata,ydata,'ro')
            plt.plot(plot_range, parabola(plot_range), 'b')
            plt.plot(plot_range, function(plot_range), 'k')
            plt.show()

        x_i = (r+s)/2.0 - (fs-fr)*(t-r)*(t-s) / 2 / ((s-r)*(ft-fs) - (fs-fr)*(t-s))

        x[k+3] = x_i
        y[k+3] = function(x[k+3])
        err = abs(y[k+3] - y[k+2])
        k += 1
    print(x)
    return x[k+2]


def line_search(function, gradient, x_k, alpha_init, beta_1, beta_2):
    # function: the objective function
    # gradient: the gradient function
    # x_k: current point of iteration
    # alpha_init: a value of alpha to start with
    # beta_1: parameter for the 1st Wolfe condition
    # beta_2: parameter for the 2nd Wolfe condition

    # Check Inputs
    if not 0 < beta_1 < beta_2 < 1:
        raise ValueError("To ensure convergence, pick 0 < beta_1 < beta_2 < 1.")

    # Initialization
    i = 0
    alpha = alpha_init
    lower = 0.
    upper = 2*alpha_init

    # 1st Wolfe condition
    gradxk = gradient(x_k)
    normgradxk = np.linalg.norm(gradxk)
    newx = x_k - np.dot(alpha,gradxk)
    wolfe_1 = function(newx) <= function(x_k) - alpha * beta_1 * normgradxk

    # 2nd Wolfe condition
    wolfe_2 = np.linalg.norm(gradient(newx),2) <= beta_2 * normgradxk

    while not (wolfe_1 and wolfe_2):

        if not wolfe_1:
            upper = alpha
            alpha = 0.5 * (lower + alpha)

        elif not wolfe_2:
            lower = alpha
            alpha = 0.5 * (alpha + upper)

        # Check both Wolfe conditions
        # 1st Wolfe condition
        newx = x_k - np.dot(alpha, gradxk)
        wolfe_1 = function(newx) <= function(x_k) - alpha * beta_1 * normgradxk
        # 2nd Wolfe condition
        wolfe_2 = np.linalg.norm(gradient(newx)) >= beta_2 * normgradxk

        """
        print("wolfe1",wolfe_1)
        print("wolfe2",wolfe_2)
        print("alpha",alpha)
        print("newx", newx)
        """

        i += 1
        if i >1000:
            raise ValueError("Failed to converge.")

    return newx


def gradient_descent(function, gradient, x_init=[], tolerance= 1.0e-6, maxIter=100, ls_method="backtracking", step_size=0.0):

    # function: the objective function, a n-variable function that takes in an array and outputs a real number
    # gradient: the gradient function of f, a list of n partial derivatives of f with respect to the x_i's
    # x_init: the initial approximation, an array that matches the dimension of f, default is the zero vector
    # ls_method: select backtracking for backtracking line search, fixed for a fixed step size
    # step_size: size of one step, required if choose fixed step

    # Check inputs
    x_init = np.asarray(x_init)

    # Initialization
    x = x_init
    gradx = gradient(x_init)
    err = tolerance * 100
    k = 0

    # Fixed step size
    if ls_method == "fixed":
        if step_size <= 0.0:
            raise ValueError("Enter a positive step size.")

        while err > tolerance and k <= maxIter:
            gradx = gradient(x)
            err = num.norm(gradx,2)
            change = step_size * gradx
            x = x - change
            k += 1

    # Inexact line search with backtracking
    elif ls_method == "backtracking":
        alpha_0 = 0.1   # I don't know how to make an educated guess here
        while err > tolerance and k <= maxIter:
            gradx = gradient(x)
            err = np.linalg.norm(gradx)
            k += 1
            x = line_search(function, gradient, x, alpha_0, beta_1=1.0e-3, beta_2=1-0.8)
            print(k,x)

    else:
        raise ValueError("Choose a valid line search method: ls for line search using SPI; fixed for a fixed step size.")

    optimality = function(x)

    return x, optimality

"""
    Exact line search with SPI (I can't get this to work yet)
    elif ls_method == "SPI":
        while err > tolerance and k <= maxIter:
            gradx = gradient(x)
            err = num.norm(gradx, 2)
            linef = lambda alpha: function(x - alpha * gradx)
            x = SPI(linef, x, tolerance=1000 * tolerance, maxIter=20)
            k += 1
"""
