MultilayerEM
------------
------------

Owen D. Miller
Oct./Nov. 2013
 
A set of c/c++ (more like c + classes) codes for a variety of electromagnetics
computations for multilayer stacks of metals/dielectrics.  The code uses a stable 
scattering matrix computation, in contrast to the more typically used (and 
numerically unstable) transmission matrix formalism.  For a stack composed of 
arbitrarily many layers, one can compute:

- Scattering matrices
- Reflection / Transmission (field and power coefficients)
- Density of states (electric DOS) at any point in the stack
- Radiative heat transfer, between either a single point or a whole layer and 
a user-specified set of flux planes.

The heat transfer is computed through the fluctuation-dissipation theorem, which 
provides a simple expression for the stochastic sources in media at a given
temperature. 

A ref. for both the scattering matrix approach and the heat transfer formulation
is:

Francoeur et. al., "Solution of near-field thermal radiation in one-dimensional
 layered media using dyadic Green's functions and the scattering matrix method,"
 Journal of Quantitative Spectroscopy and Radiative Transfer 110, 2002 (2009)
[doi](http://dx.doi.org/10.1016/j.jqsrt.2009.05.010)

The codes are written in standard non-dimensional units.  Simple conversion 
factors are commented around the corresponding functions.

Codes:
* heat\_transfer\_sim.cpp: a working example of how to compute emission rates 
in a large complex multilayer system.

* dos\_sim.cpp: similarly, a working example of how to compute the density of 
states.

* SMatrix.cpp: a class to hold scattering matrices of a structure at a given 
frequency (w) and parallel wavevector (kp)

* mlgeo.cpp: a multilayer geometry class, encapsulating the permittivity and
thickness of each layer (including the embedding half-spaces)

* fieldFunctions.cpp: contains the functions to compute flux rates, the density 
of states, and the reflection/transmission coefficients.  Also contains 
auxiliary functions to compute blackbody flux rates and the density of states 
in vacuum (in non-dimensional units, again).

* materials.cpp: (NOT REQUIRED) Note that this is an auxiliary code I have 
provided.  It is a code I use with functions for common material models
(in the infrared).  Users can specify permittivities however they see fit.

* numericalIntegration.cpp: (NOT REQUIRED) This file provides functions to 
perform numerical integrations over wavevector and possible frequency.  
It uses the [Cubature](ab-initio.mit.edu/wiki/index.php/Cubature) package
written by [Steven G. Johnson](math.mit.edu/~stevenj)

* Many of the above have corresponding header files

Note that as currently written, heat\_transfer\_sim.cpp and dos\_sim.cpp do 
each use a lambda, so I compile them with the -std=c++11 flag.  However, this 
is not at all required (it is only to comply with my own material function 
definitions) and no other aspect of the code depends on c++11.

On Windows, you may have to compile with the flag -D_USE_MATH_DEFINES
