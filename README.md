# Modal Discontinuous Galerkin Method in MATLAB

This repository contains the modal discontinuous Galerkin method implemented in MATLAB/Octave. Note that the purpose of this code is to test specific setup and methods so the code is not optimized nor fast and the emphasis is mainly on its readability and simplicity. 

Folder `./examples` contains a few scripts putting the building blocks together.

### Done so far:
- Gauss-Legendre and Lobatto-Gauss-Legendre quadrature rules
- Legendre basis functions
- 1D implementation
- Periodic boundary conditions only
- Linear scalar advection equation and Burgers' equation
- Explicit time-integration (Strong Stability-Preserving Runge-Kutta methods)

### TODO:
- Implicit time-integration
- Steady problems
- Hybridized DG
- Diffusion discretization
- Euler equations
- Stabilization of discontinuous solutions