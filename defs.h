/*		A S Sharma June 2011		*/


/*				DEFINITIONS REQUIRED BY distribution.h				*/
// DIMS is dimension of problem / velocities
#define DIMS 2
// EPS is threshold to truncate distribution
#define EPS 1e-6
// RHS parameters for Lorentz
#define B 1
#define SIGMA 4
#define R 48
// diffusion variance when stepping distribution with diffusion enabled
#define VAR 0

// for the fixed grid, this is the grid size in boxes
#define GRIDSIZE 400
