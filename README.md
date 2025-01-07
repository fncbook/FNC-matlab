# FNC-matlab

These are the MATLAB versions of the functions used in the book [*Fundamentals of Numerical Computation*](https://fncbook.github.io/fnc) by Tobin A. Driscoll and Richard J. Braun. There have been some implementation changes since the print edition of the text.

- `hatfun` now returns a callable function of $x$ rather than requiring the evaluation points at the time of construction.
- The IVP solvers use MATLAB's [new IVP structure](https://www.mathworks.com/help/releases/R2024b/matlab/ref/ode.html) to set up the problem. This makes it potentially less confusing to pass parameters into the ODE function.
- The boundary conditions of a BVP in Chapter 10 are specified a little differently.
- There is a `tensorgrid` function that creates several values and functions useful to the discussion and algorithms of Chapter 13.
- The `newtonpde` function has been replaced by a more general `elliptic` function to solve an elliptic PDE on a rectangular domain.

## Installation

Download the latest release as a zip file and unpack somewhere on your computer. Add the `FNC-matlab` directory to your MATLAB path. By typing `addpath /path/to/FNC-matlab` at the MATLAB prompt, or by typing `pathtool` and using the GUI. The `pathtool` method allows you to save the path for all future MATLAB sessions.