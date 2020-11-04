# Numerical Analysis by Mingxuan He
## Description:
I'm currently (Fall 2020 term 2) taking MAT-317 Numerical Analysis, an advanced-level applied mathematics course at Grinnell College. During the course I would regularly implement numerical algorithms in Python and (in some rare occasions) in MATLAB, and I will be sharing my code on Github along the way. 

Environment: Python 3.7 with Jetbrains Pycharm

Packages used: numpy, matplotlib

Citations: Some code snippets and most test files are taken from my instructor prof. Jeffery Blanchard.

## Algorithms:
Below is a brief description of what each of the algorithms in the package does. I did not have time for a comprehensive documentation, but you can refer to comments in the code for more information.

### Nested Polynomials (*polynest*)
- Implemented in numerics0_he
- This algorithm evaluate a polynomial with real coefficients a and basepoints b at the point x and returns that value as y
- Converting a polynomial into the nested form before evaluating is much more efficient than regular methods since it requires significantly less operations.

### Bisection method for root finding (*bisection*)
- Implemented in numerics1_he
- This algorithm finds a root of a function by bisecting a real interval which we know contains the root.

### Fixed point iteration (*fixedpt*)
- Implemented in numerics1_he
- This algorithm finds a fixed point of a function by iteratively taking its own result as the function input
- Requires the function and initial value to satisfy certain conditions in order to give a convergent fixed point (see Thm 1.2.a and Thm 1.2.b)
