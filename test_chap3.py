import numpy as np
import matplotlib.pyplot as plt
import numerics3_he as num


divided_diff = 0
interp_func = 0
run_plotting = 0
run_cheby = 1
run_cubicspline = 0  # Not working just yet.

# The plotting depends on the function we create.
if (run_plotting):
    interp_func = 1

if (divided_diff):
    # The data for example 3.6 in the text.
    data = np.array([[-1.0, 0, 2, 3],[-5.0, -1,  1, 11]])
    d = num.newtondd(data[0,:],data[1,:])
    print(d)

if (interp_func):
    # A problem taken from Burden and Faires
    x = [1.0, 1.3, 1.6, 1.9, 2.2]
    y = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]

    fxn = num.newtonInterp(x,y)
    print("f(1.5)=%f, f(2.0) = %f" % (fxn(1.5),fxn(2.0)))

if (run_plotting):
    t = np.arange(-5,5,0.001)

    # plotted the function
    plt.plot(t, fxn(t), 'r')
    # plot the x-axis
    plt.plot(t, 0.0*t, 'k')
    plt.plot(x,y,'bo')
    plt.show()

    t = np.arange(0.5,2.5,0.01)

    # plotted the function
    plt.plot(t, fxn(t), 'r')
    # plot the x-axis
    plt.plot(t, 0.0*t, 'k')
    plt.plot(x,y,'bo')
    plt.show()

if (run_cheby):
    numrts = 10
    rts = num.chebyshevRoots(numrts)
    print("roots",rts)
    intvl = [0.0,np.pi]
    print("range",range)

    func = lambda x: np.sin(x)
    chebfunc = num.chebyshevInterp(func, numrts, intvl)

    def mysin(x):
        s = 1
        temp = np.mod(x,2*np.pi)
        if (temp>np.pi):
            temp = 2*np.pi-temp
            s = -1
        elif (temp > np.pi/2.0):
            temp = np.pi - temp
        return s*chebfunc(temp)

    examp_pts = np.array([.01,1.0, 2.0, 3.0, 4.0])
    print("The input followed by the values of numpy's sine function, our version, and the chebyshev interpolating polynomial.")
    for idx in range(len(examp_pts)):
        xval = examp_pts[idx]
        print("%d\t %0.5f\t %0.5f\t  %0.5f" % (xval,func(xval),mysin(xval),chebfunc(xval)))

    t = np.arange(intvl[0], intvl[1], 0.001)

    # plotted the function
    plt.plot(t, np.sin(t), 'r')
    plt.plot(t,chebfunc(t),"y")
    # plot the x-axis
    plt.plot(t, 0.0 * t, 'k')
    plt.plot(examp_pts,chebfunc(examp_pts), 'bo')
    plt.show()


if (run_cubicspline):
    x = [0,1,2]
    y = [3,-2,1]
    print("cubes")
    ccc = num.cubiccoeff(x,y)
    print(ccc)

