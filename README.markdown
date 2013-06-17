Grid-based Estimation, Exploiting Sparsity
==========================================
*GBEES v1.1.4*
--------------

Maintainer
----------

A S Sharma  
Aerodynamics and Flight Mechanics  
Faculty of Engineering and the Environment  
University of Southampton, UK  
http://www.personal.soton.ac.uk/as5v12/


Citation
--------

Efficient grid-based Bayesian estimation of nonlinear low-dimensional systems with sparse non-Gaussian PDFs  
T R Bewley, A S Sharma  
Automatica 48 (7), 2012, pp. 1286-1290, [DOI](http://dx.doi.org/10.1016/j.automatica.2012.02.039)  

Also available on [ArXiV](http://arxiv.org/abs/1301.4866v1) and the maintainer's [personal web site](http://www.personal.soton.ac.uk/as5v12/).


Abstract
--------

Bayesian estimation strategies represent the most fundamental formulation of the state estimation problem available, and apply readily to nonlinear systems with non-Gaussian uncertainties. The present paper introduces a novel method for implementing grid-based Bayesian estimation which largely sidesteps the severe computational expense that has prevented the widespread use of such methods. The method represents the evolution of the probability density function (PDF) in phase space, $p_{\x}(\x',t)$, discretized on a fixed Cartesian grid over _all_ of phase space, and consists of two main steps: (i) Between measurement times, $p_{\x}(\x',t)$ is evolved via numerical discretization of the Kolmogorov forward equation, using a Godunov method with second-order corner transport upwind correction and a total variation diminishing flux limiter; (ii) at measurement times, $p_{\x}(\x',t)$ is updated via Bayes' theorem. Computational economy is achieved by exploiting the localised nature of $p_{\x}(\x',t)$. An ordered list of cells with non-negligible probability, as well as their immediate neighbours, is created and updated, and the PDF evolution is tracked *only* on these active cells.


The code
--------

The code is in C++ and spits out files for munging by Matlab post-processing. It should compile with any recent version of gcc or llvm without too much complaining. Example outputs may be seen in the paper.

The code can handle N-dimensional systems (not just 3D) however no guarantees for performance are given for isystems of very high dimension.


Usage
-----

1. hack `defs.h` to reflect any needed parameters
2. hack `distribution::u()` in GBEES.cpp to reflect your function
3. compile & run

You can see examples in the code for solid body rotation, and for the Lorenz equations case from the paper.


Changelog
---------

* v0.1 used vectors; slow
* v0.2 uses maps; O(log(n)) most ops
* v0.3 implements key by {i,j,k} and sparserotate neighbour tracking by iterators, uses lists efficiently, Gudenov, flux limiting etc
* v0.3.1 CFL targeting & bugfixes
* v1.0.1 bugfixes, truth model
* v1.0.2 bugfixes
* v1.0.3 start of noise
* v1.0.4 Gaussian blur implemented
* v1.0.6 Entropy measures - a first attempt
* v1.0.7 Minor speed improvements; bug fix
* v1.0.8 Bayesian measurement updates; if statements rewritten
* v1.0.9 Gaussian reference distribution for m im H calculation; and dump/plot distrib every given time
* v1.0.10 Rewrote `distribution::drop()` to not loop over all elements -- 35% speed gain!
* v1.0.11 Changed `update()` to unchecked insert, and consolidate
* v1.1.0 Now uses `std::map`; this gives O(n log(n)) performance -- *MUCH* faster -- now critical step is timestepping
* v1.1.1 fixed bug in second order flux correction terms
* v1.1.2 things properly arranged in header files and .cpp files
* v1.1.3 coded 2d leveque test case
* v1.1.4 put up on GitHub


Future features
---------------

- Bugfix so that singular iterators are not copied
- More elegantly implement dynamics function `distribution::u()`
- Better output for plotting; real-time hooks into paraview instead?
- Parallelisation / OpenMP? -- may not be parallelisable
- Account for P truncation error
- Particle filtering
- Adaptive gridding (instead of sparsity)


Licensing
---------

The code is available under a modified MIT licence. See separate LICENSE file for details.
