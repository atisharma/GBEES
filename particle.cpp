/*		A S Sharma June 2009		*/

#include <particle.h>

/**********************************************/
//	particle CLASS
/**********************************************/

particle::particle(string fname, distribution *P, double x0[]) {
	for (int n=0;n<DIMS;n++) state[n]=x0[n];
	particle_file.open(fname.c_str(), ios::out);
	dist = P;
}

particle::~particle() {
	particle_file.close();
}

void particle::step(double dt) {
	double temp[DIMS];
	double k1[DIMS]; double k2[DIMS];
	double k3[DIMS]; double k4[DIMS];
	
	for (int n=0;n<DIMS;n++) k1[n]=(dist->u)(state,n);
	for (int n=0;n<DIMS;n++) temp[n]=state[n]+k1[n]*dt/2;
	for (int n=0;n<DIMS;n++) k2[n]=(dist->u)(temp,n);
	for (int n=0;n<DIMS;n++) temp[n]=state[n]+k2[n]*dt/2;
	for (int n=0;n<DIMS;n++) k3[n]=(dist->u)(temp,n);
	for (int n=0;n<DIMS;n++) temp[n]=state[n]+k3[n]*dt;
	for (int n=0;n<DIMS;n++) k4[n]=(dist->u)(temp,n);
	for (int n=0;n<DIMS;n++) state[n] += dt*(k1[n]/6 + (k2[n] + k3[n])/3 + k4[n]/6);
}

void particle::save(double t) {
	for (int n=0;n<DIMS;n++){
		particle_file << state[n] << "\t";
	}
	particle_file << t;
	particle_file << endl;
}
