/*		A S Sharma June 2009		*/

#ifndef PARTICLE_H
#define PARTICLE_H


#include <defs.h>
#include <distribution.h>

/**********************************************/
//	particle CLASS
/**********************************************/
class particle {
private:
	ofstream particle_file;
	distribution * dist;
public:
	double state[DIMS];
	particle(string, distribution *, double[]);
	~particle();
	void step(double dt);
	void save(double t);
};

#endif
