import numpy as np
import numerics_he as num
import numericsproj_he as proj
import matplotlib.pyplot as plt

test_SPI = 0
test_GD_unconstrained = 0
test_GD_constrained = 1

if test_SPI == 1:

    print("Testing for SPI:")

    def testf1(x):
        return x**6. - 11*x**3 + 17*x**2 - 7*x + 1

    x_inits = [0,2,10]

    print("Solution:")
    print(proj.SPI(testf1, x_inits, maxIter=15, plot=True, plot_range=np.arange(0,1,0.01)))


if test_GD_unconstrained == 1:

    print("Testing for Gradient Descent (unconstrained):")

    def testf2(xy):
        x = xy[0]
        y = xy[1]
        return 5.*x**4 + 4*x**2*y - x*y**3 + 4*y**4 - x

    def testf2_grad(xy):
        x = xy[0]
        y = xy[1]
        return np.array([20.*x**3 + 8*x*y - y**3 - 1, 4*x**2 - 3*x*y**2 + 16*y**3])


    x_inits = np.array([1.,-1.])


    step_sizes = [0.01, 0.05, 0.1]
    for ss in step_sizes:
        x, optimality = proj.gradient_descent(testf2, testf2_grad, x_inits, ls_method="fixed", step_size=ss)
        print("With fixed step size:", ss)
        print("x: ", x)
        print("min f(x): ", optimality)

    x, optimality = proj.gradient_descent(testf2, testf2_grad, x_inits, ls_method="backtracking")

    #x,optimality = proj.gradient_descent(mysqr,derquad,np.array([3.]),ls_method="backtracking")

    print("Backtracking results:")
    print("x: ", x)
    print("min f(x): ", optimality)

    # testing line search only
    #x = proj.line_search(testf2, testf2_grad, x_inits, 15., beta_1=1.0e-3, beta_2=0.8)

    def mysqr(x):
        return x**2

    def derquad(x):
        return 2*x

    #x = proj.line_search(mysqr, derquad, np.array([3.]), 1., beta_1=1.0e-3, beta_2=0.8)



if test_GD_constrained == 1:

    # the objective
    def testf3(xy):
        x = xy[0]
        y = xy[1]
        return -x**(1/3)*y**(2/3)

    # the constraint
    def testg1(xy):
        x = xy[0]
        y = xy[1]
        return x+2*y-10

    # the lagrangian
    def testf3_lag(xylam):
        x = xylam[0]
        y = xylam[1]
        lam = xylam[2]
        return testf3(np.array([x,y])) - lam * testg1(np.array([x,y]))

    # gradient of the lagrangian
    def testf3_lag_grad(xylam):
        x = xylam[0]
        y = xylam[1]
        lam = xylam[2]
        return np.array([-(1/3)*x**(-2/3)*y**(2/3)-lam,
                         -(2/3)*x**(1/3)*y**(-1/3)-lam,
                         -testg1(np.array([x, y]))])

    #print("Constrained optimum with Lagrangian:")
    #x_inits = np.array([3, 3., 1])
    #x, optimality = proj.gradient_descent(testf3_lag, testf3_lag_grad, x_inits, ls_method="fixed", step_size=0.1)
    #x, optimality = proj.gradient_descent(testf3_lag, testf3_lag_grad, x_inits, ls_method="backtracking", tolerance=1.0e-3, maxIter=100)
    #print("x",x)
    #print("min f(x)",optimality)

    # the obj with penalty
    def testf3_pen(xy):
        x = xy[0]
        y = xy[1]
        return testf3(xy) + 100*testg1(xy)**2

    # gradient of the obj with penalty
    def testf3_pen_grad(xy):
        x = xy[0]
        y = xy[1]
        return np.array([0-(1/3) * x ** (-2/3) * y ** (2/3) + 100*(2*x+4*y-20),
                         0-(2/3) * x ** (1/3) * y ** (-1/3) + 100*(4*x+8*y-40)])

    print("Constrained optimum with penalty:")
    x_inits = [5.,5.]
    x, optimality = proj.gradient_descent(testf3_pen, testf3_pen_grad, x_inits, ls_method="backtracking", tolerance=1.0e-3,maxIter=100000)
    print("x",x)
    print("min f(x)",optimality)

"""    # the objective
    def testf2(xy):
        x = xy[0]
        y = xy[1]
        return 5.*x**4 + 4*x**2*y - x*y**3 + 4*y**4 - x

    # the constraint
    def testg2(xy):
        x = xy[0]
        y = xy[1]
        return x+y-1

    # the lagrangian
    def testf2_lagrangian(xylambda):
        x = xylambda[0]
        y = xylambda[1]
        lam = xylambda[2]
        return testf2([x,y]) - lam * testg2([x,y])

    def testf2_lag_grad(xylambda):
        x = xylambda[0]
        y = xylambda[1]
        lam = xylambda[2]
        return np.array([20.*x**3 + 8*x*y - y**3 - 1 - lam, 4*x**2 - 3*x*y**2 + 16*y**3 - lam, x + y-1])

    x_inits = np.array([1.,-1., 1])

    x, optimality = proj.gradient_descent(testf2_lagrangian, testf2_lag_grad, x_inits, ls_method="backtracking")

    print("Constrained optimum with Lagrangian:")
    print("x",x)
    print("min f(x)",optimality)


    def testf2_penalty(xy):
        x = xy[0]
        y = xy[1]
        return testf2(xy) + 10*testg2(xy)**2

    def testf2_penalty_grad(xy):
        x = xy[0]
        y = xy[1]
        return np.array([20.*x**3 + 8*x*y - y**3 - 1 +20*x+20*y-20, 4*x**2 - 3*x*y**2 + 16*y**3 +20*x+20*y-2])

    x_inits = [1,-1.]

    x, optimality = proj.gradient_descent(testf2_penalty, testf2_penalty_grad, x_inits, ls_method="backtracking")

    print("Constrained optimum with penalty:")
    print("x",x)
    print("min f(x)",optimality)
"""
