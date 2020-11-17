# Numerical Analysis by Mingxuan He
## Description:
I'm currently (Fall 2020 term 2) taking MAT-317 Numerical Analysis, an advanced-level applied mathematics course at Grinnell College. During the course I would regularly implement numerical algorithms in Python (maybe in MATLAB in some rare occasions), and I will be sharing my code on Github along the way. 

Environment: Python 3.7 with Jetbrains Pycharm

Packages used: numpy, matplotlib

Citations: Some code snippets and most test files are taken from my instructor prof. Jeffery Blanchard.

## Algorithms:
Below is a brief description of what each algorithm in the package does. I did not have time for a comprehensive documentation, but you can refer to comments in the code for more information.

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

### LU factorization (*lu, lusolve, noswaplu, noswaplusolve*)
- Implemented in numerics2_he
- Usage: Precise root finding for systems of linear equations 
- This algorithm finds the root for a system of linear equations (in matrix form, Ax=b where A is a n by n matrix and b is a vector of size n) by factoring the permuted coefficient matrix PA (alternatively A itself) into a lower-triangular matrix L and an upper triangular matrix U, then solve PAx=LUx=Pb (alternatively Ax=LUx=b) using forward-backward substitution.
- Helper functions: 
  - Row operations: *rowswap, rowscale, rowdiff*
  - Foward and backward substitution: *forwardsub, backwardsub, fbsolve*

### Jacobi Method (*jacobi*)
- Implemented in numerics2_he
- Usage: Iterative rooting finding for systems of linear equations
- This algorithm solves Ax=b iteratively using Jacobi Method (additive decomposition)
- Helper function: *add_decomp*

### Gauss-Siedel Method (*gausssiedel*)
- Implemented in numerics2_he
- Usage: Iterative rooting finding for systems of linear equations
- This algorithm solves Ax=b iteratively using Gauss-Siedel Method (additive decomposition)

### Least Squares (*leastSquares_lu*)
- Implemented in numerics4_he
- Usage: Modeling data
- This algorithm uses LU factorization to solve the system A^tAx=A^tb to find a best-fit model that satisfies the least squares condition.
