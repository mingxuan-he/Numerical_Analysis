import numerics0_he as num0
import numerics1_he as num1
import numerics2_he as num2
import numerics3_he as num3
import numerics4_he as num4
import numerics5_he as num5


def polynest(x, a, b=[]):
    return num0.polynest(x, a, b=[])


def bisection(function, left_endpoint, right_endpoint, tolerance=1.0e-6, maxIter=100):
    return num1.bisection(function, left_endpoint, right_endpoint, tolerance=1.0e-6, maxIter=100)


def fixedpt(function, xinit, tolerance=1.0e-6, maxIter=100):
    return num1.fixedpt(function, xinit, tolerance=1.0e-6, maxIter=100)


def newton(function, dfunction, xinit, tolerance=1.0e-6, maxIter=100):
    return num1.newton(function, dfunction, xinit, tolerance=1.0e-6, maxIter=100)


def secant(function, xinit, tolerance=1.0e-6, maxIter=100):
    return num1.secant(function, xinit, tolerance=1.0e-6, maxIter=100)


def norm(vec, l):
    return num2.norm(vec, l)


def rowdiff(matrix,k,j,scale=1.0):
    return num2.rowdiff(matrix,k,j,scale=1.0)


def rowswap(matrix,k,j):
    return num2.rowswap(matrix,k,j)


def rowscale(matrix,k,scale=1.0):
    return num2.rowscale(matrix,k,scale=1.0)


def noswapLU(matrix):
    return num2.noswapLU(matrix)


def backsub(U,b):
    return num2.backsub(U,b)


def forwardsub(L,b):
    return num2.forwardsub(L,b)


def fbsolve(L,U,b):
    return num2.fbsolve(L,U,b)


def noswaplusolve(A,b):
    return num2.noswaplusolve(A,b)


def lu(matrix):
    return num2.lu(matrix)


def lusolve(A,b):
    return num2.lusolve(A,b)


def add_decomp(A):
    return num2.add_decomp(A)


def jacobi(A, b, xinit, tolerance=1.0e-6, maxIter=100):
    return num2.jacobi(A, b, xinit, tolerance=1.0e-6, maxIter=100)


def gausssiedel(A, b, xinit, tolerance=1.0e-6, maxIter=100):
    return num2.gausssiedel(A, b, xinit, tolerance=1.0e-6, maxIter=100)


def newtondd(xdata,ydata):
    return num3.newtondd(xdata,ydata)


def newtonInterp(xdata,ydata):
    return num3.newtonInterp(xdata,ydata)


def cubiccoeff(xdata, ydata, end_condition="natural", df = lambda x: 0):
    return num3.cubiccoeff(xdata, ydata, end_condition="natural", df = lambda x: 0)


def cubic_eval(xdata, ydata, x, end_condition="natural", df = lambda x: 0):
    return num3.cubic_eval(xdata, ydata, x, end_condition="natural", df = lambda x: 0)


def chebyshevRoots(num_roots):
    return num3.chebyshevRoots(num_roots)


def chebyshevInterp(function, num_roots, interval=[-1,1], interp_method="Newton"):
    return num3.chebyshevInterp(function, num_roots, interval=[-1,1], interp_method="Newton")


def leastSquares_lu(A,b):
    return num4.leastSquares_lu(A,b)


def householder(vec):
    return num4.householder(vec)


def qr(matrix):
    return num4.qr(matrix)


def qrsolve(A,b):
    return num4.qrsolve(A,b)


def simpson(function, a, b):
    return num5.simpson(function, a, b)


def compositeSimpson(function, n, a, b):
    return num5.compositeSimpson(function, n, a, b)


def adaptiveSimpson(function, a, b, tolerance=1.0e-12, Sab=0.0, recursion_counter=0, max_depth=20):
    return num5.adaptiveSimpson(function, a, b, tolerance=1.0e-12, Sab=0.0, recursion_counter=0, max_depth=20)


def gaussQuad(function,n,a,b):
    return  num5.gaussQuad(function,n,a,b)



