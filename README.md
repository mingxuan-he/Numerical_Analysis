# Numerical Algorithms Package by Mingxuan He
## Description:
I am a mathematics and economics double major at Grinnell College. Currently (Fall 2020 Term 2) I'm taking MAT-317 Numerical Analysis, an advanced-level applied mathematics course. Throughout the course I will regularly implement numerical algorithms in Python (maybe MATLAB in some rare occasions), and I will be sharing my code on Github along the way. Feel free to use this package however you'd like as long as this github repo is cited.

Environment: Python 3.7 with Jetbrains Pycharm Ver. 2020.1

Packages used: numpy, matplotlib

Citations: Some code snippets and most test files are taken from my instructor prof. Jeffery Blanchard.

### Personal Info:
Linkedin: https://www.linkedin.com/in/mingxuanhe/

Personal Website: https://mingxuanhe.godaddysites.com/

Feel free to reach out for questions/comments/concerns/suggestions!

## Algorithms:
Below is a brief description of what each algorithm in the package does. I did not have time for a comprehensive documentation, but in most cases the comments in the code will provide more detailed information.

### Nested Polynomials (*polynest*)
- Implemented in numerics0_he
- Usage: Evaluating polynomials
- This algorithm evaluate a polynomial with real coefficients a and basepoints b at the point x and returns that value as y
- Converting a polynomial into the nested form before evaluating is much more efficient than regular methods since it requires significantly less operations.

### Bisection method (*bisection*)
- Implemented in numerics1_he
- Usage: Root finding
- This algorithm finds a root of a function by bisecting a real interval which we know contains a root.
- Requirement: an interval (a,b) that is known to contain a root.

### Fixed point iteration (*fixedpt*)
- Implemented in numerics1_he
- Usage: Finding a fixed point x such that f(x)=x
- This algorithm finds a fixed point of a function by iteratively taking its own result as the function input
- Requirement: the function and initial value likely but not necessarily satisfy certain conditions in Thm 1.2.a and Thm 1.2.b

### Newton's Method (*newton*)
- Implemented in numerics1_he
- Usage: Root finding
- This algorithm finds a root of a function by iteratively linearizing the function (using its derivative) to approach the root.
- Requirement: the 1st derivative of the function; r is likely but not necessarily a simple root.

### Secant Method (*secant*)
- Implemented in numerics1_he
- Usage: Root finding (derivative-free)
- This algorithm finds a root of a function by iteratively linearizing the function (using its secant) to approach the root. It's similar to Newton's Method, except that the secant method only requires evaluating f(x), but not its derivative. Less accurate but useful when the derivative is unknown.
- Requirement: r is likely but not necessarily a simple root.

### LU factorization (*lu, lusolve*)
- Implemented in numerics2_he
- Usage: Precise root finding for systems of linear equations 
- This algorithm finds the root for a system of linear equations (in matrix form, Ax=b where A is a n by n matrix and b is a vector of size n) by factoring the permuted matrix PA (permutation is done using partial pivoting) into a lower-triangular matrix L and an upper triangular matrix U, then solve PAx=LUx=Pb using forward-backward substitution.
- Helper functions: 
  - A=LU factorization (no swap): *noswaplu, noswaplusolve*
  - Row operations: *rowswap, rowscale, rowdiff*
  - Foward and backward substitution: *forwardsub, backwardsub, fbsolve*

### Jacobi Method (*jacobi*)
- Implemented in numerics2_he
- Usage: Iterative rooting finding for systems of linear equations
- This algorithm solves Ax=b iteratively using Jacobi Method (additive decomposition)
- Helper functions: *add_decomp*, *norm*
- The algorithm will converge to a unique root if A is strictly diagonally dominant.

### Gauss-Siedel Method (*gausssiedel*)
- Implemented in numerics2_he
- Usage: Iterative rooting finding for systems of linear equations
- This algorithm solves Ax=b iteratively using Gauss-Siedel Method (additive decomposition)

### Least Squares (*leastSquares_lu*)
- Implemented in numerics4_he
- Usage: Modeling data (approximation)
- This algorithm uses LU factorization to solve the system A^tAx=A^tb to find a best-fit model that satisfies the least squares condition.

### QR Factorization with Householder Reflectors (*qr, qrsolve*)
- Implemented in numerics4_he
- Usage: Solving systems of equations (can be used for least squares)
- This algorithm solves Ax=b where A is an m by n matrix. Uses householder reflectors to factor A into an orthogonal matrix Q and an upper triangular matrix R, then solves Rx=Q^Tb
- Helper function: *householder*

### Newton Polynomial Interpolation (*newtonInterp*)
- Implemented in numerics3_he
- Usage: Modeling data (precise interpolation)
- This algorithm models data points (x_i,y_i) with a Newton Polynomial using a divided differences table.
- Helper function: *newtondd*

### Chebyshev Interpolation (*chebyshevRoots,chebyshevInterp*)
- Implemented in numerics3_he
- Usage: Approximate a function on an interval
- This algorithm approximates a function on a given interval by sampling on Chebyshev nodes (roots of the nth Chebyshev polynomial), thus minimizing interpolation errors.
