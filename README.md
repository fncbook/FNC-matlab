# FNC-matlab

These are the MATLAB versions of the functions used in the book [*Fundamentals of Numerical Computation*](https://fncbook.github.io/fnc) by Tobin A. Driscoll and Richard J. Braun. There have been some implementation changes since the print edition of the text.

- `hatfun` now returns a callable function of $x$ rather than requiring the evaluation points at the time of construction.
- The boundary conditions of a BVP in Chapter 10 are specified a little differently.
- There is a `tensorgrid` function that creates several values and functions useful to the discussion and algorithms of Chapter 13.
- The `newtonpde` function has been replaced by a more general `elliptic` function to solve an elliptic PDE on a rectangular domain.

## Installation

- To install in MATLAB Online, you can click this button:
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=fncbook/FNC-matlab&file=toolbox/release/FNC.mltbx)
- To install the functions locally, you can download the `mltbx` file from the [releases page](https://github.com/fncbook/FNC-matlab/releases). Double-click the file, or load it in to the MATLAB Add-On Manager.
- The old-school way is to download a source code zip file from the [releases page](https://github.com/fncbook/FNC-matlab/releases) and unpack it somewhere on your computer. Then add its top-level directory to your MATLAB path by using the `addpath` command or the `pathtool` GUI.
