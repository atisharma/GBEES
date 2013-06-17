/*		A S Sharma 6-10/2009				*/
/*		originally based on dbees.m         */
/*      by T Bewley 6/2009	            	*/
//

/*				DEFINITIONS REQUIRED BY CODE IN defs.h			*/
#include <defs.h>

#include <iostream>
#include <sstream>
#include <distribution.h>
#include <particle.h>
#include <ctime>

using namespace std;

//-------------------------------------------- u - defines RHS of particle
//		upper case: Lorenz; lower case: solid body rotation
double distribution::u(double x[DIMS], int dim) {
	switch (dim) {
		case 0:
			//return SIGMA*(x[1]-x[0]);
			return x[1];
		case 1:
			//return -x[1] -x[0]*x[2];
			return -x[0];
		//case 2:
			//return -B*x[2] + x[0]*x[1] - B*R;
			//return -0.5;
	}
	//	throw an exception
	assert (DIMS==2);
	return 0;
}


//-------------------------------------------- randn()
// Marsaglia polar method of generating Gaussian random function
double randn(double variance) {
	double r=100.0, x1, x2;
	while (!(r<1)) {
		x1=((double)rand()/(double)RAND_MAX)*2 - 1; // [-1,1]
		x2=((double)rand()/(double)RAND_MAX)*2 - 1; // [-1,1]
		r = x1*x1 + x2*x2;
	}
	return x2*sqrt(-2*variance*log(r)/r);
}

double randn() {
	return randn(1);
}

/**********************************************/
int main() {

	double const dx=0.01; // LeVeque test case
	//double const dx=0.25; //Lorenz case, between .1 and 1 is OK - watch EPS!
	double dt_max=0.001;
	double dt=.0005;
	double T=1;
	double CFL;
	double CFL_target=0.2;
	double maxv=1;
	int i=0;
	double t=0.0;
	double timer;
	double x0[3]={-11.5, -10, 9.5};
	double y0[3]={-11.5, -10, 9.5};
	string truth_file = "truth.asc";
	string frame_file = "frame.asc";
	srand(time(0)+clock());
	
	dt_max=	min((dx*dx/9)/(1e-5+VAR),dt_max);	//	auto managed by CLF; this is 3 std dev for diffusion
	
	cout << "dx=" << dx;
	cout << ", EPS=" << EPS;
	cout << ", T=" << T;
	cout << ", CFL_target=" << CFL_target;
	cout << ", dt_max=" << dt_max;
	cout << ", VAR=" << VAR;
	cout << endl;

	
	
	// LeVeque's test case ex 20.1
	// Requires DIMS=2. Ignores EPS since running different cases.
	cout << "initialising distributions" << endl;
	distribution P0(dx);
	distribution P_eps(dx);
	distribution P(dx);
	P0.test_case(); // initialise the test case
	P.test_case();
	P_eps.test_case();
	
	P0.save((string)"LeVeque_test_case.asc", 0);	//	save prob distribution
	cout << "Running GBEES..." << endl;
	timer=clock(); i=0;
	dt=0.001;
	for (t=0.0; t<2*M_PI; t+=dt) {
		i++;
		P.step(dt);
		P_eps.step(dt);
		// add new elements to distribution where appropriate
		P.update(1e-16); // very high tolerance for truncation
		P_eps.update(1e-3); // aggressive truncation
		
		// dump results out every 0.1
		if (fmod(t,0.1)<dt) {
			cout << "t=" << t;
			cout << ", P_eps size=" << P_eps.size();
			cout << ", P size=" << P.size();
			cout << ", dt=" << dt;
			maxv=P.max_v();
			CFL=maxv*dt/dx;
			cout << ", CFL=" << CFL;
			cout << endl;
		}
	}
	cout << "t=" << t << endl;
	P.save((string)"final_pure_LeVeque_test_case.asc", t);		//	save prob distribution
	P_eps.save((string)"final_eps_LeVeque_test_case.asc", t);		//	save prob distribution

	// Work out K-L divergence of normalised P0 -> P
	cout << "Distribution divergences in bits:=" << endl;
	cout << "D_KL(P0,P_eps)=" << P0.D_KL(P_eps) << endl;
	cout << "D_KL(P0,P)=" << P0.D_KL(P) << endl;
	cout << "D_KL(P0,P0)=" << P0.D_KL(P0) << endl;
	
	
	
	/*
	// The Lorenz case. Requires DIMS=3 and distribution.u() changed accordingly.
	cout << "initialising distribution" << endl;
	distribution P(dx);

	P.initialise(x0,.25);	//	Gaussian, centred on x0, variance second argument
	
	particle truth(truth_file, &P,x0);
	particle frame(frame_file, &P,y0);
	
	cout << "Sketching attractor..." << endl;
	while (t<60.0) {
		frame.step(dt);
		t+=dt; i++;
		if (!(i%10)&&(t>.2)) {
			frame.save(t);
			//	truth
		}
	}

	cout << "Running particle swarm..." << endl;
	for (int n=1;n<=200;n++) {
		t=0.0;i=0;
		stringstream swarm_name;
		swarm_name << "swarm_" << setprecision(2) << n << ".asc";
		double c0[DIMS];
		for (int d=0; d<DIMS; d++) c0[d] = x0[d] + randn(0.2);
		particle bee(swarm_name.str(), &P, c0);
		while (t<T+dt) {
			bee.save(t);
			bee.step(dt);
			t+=dt; i++;
		}
	}
	
	
	cout << "Running GBEES..." << endl;
	timer=clock(); i=0;
	for (t=0.0; t<T+dt; t+=dt) {
		i++;

		truth.step(dt);
		P.step(dt);
		
		// Bayesian update on a particular dimension (in this case 2)
		// with variance 25.0
		P.observe(truth.state[2], 25.0, 2);
		
		// add new elements to distribution where appropriate
		P.update(EPS);
		
		if (!(i%3)) {
			// check CFL
			maxv=P.max_v();
			CFL=maxv*dt/dx;
			dt=min(CFL_target*dx/maxv,dt_max);
			truth.save(t);		//	truth
		}

		// dump results out every 0.5
		if (fmod(t,0.2)<dt) {
			cout << "t=" << t;
			cout << ", size=" << P.size();
			//cout << ", timesteps=" << i;
			cout << ", dt=" << dt;
			cout << endl;
			P.save((string)"distribution.asc", t);	//	save prob distribution
		}
	}

	P.save((string)"distribution_pre_meas.asc", t);		//	save prob distribution
	P.observe(truth.state[2], 25.0, 2);	
	P.save((string)"distribution_post_meas.asc", t);	//	save prob distribution
	*/
	
	cout << "Time elapsed: " << ((double)clock()-timer)/CLOCKS_PER_SEC << "s" << endl;
	cout << "complete." << endl;
	
	return 0;
}
