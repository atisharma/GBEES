/*		A S Sharma June 2009, 2011	*/

#include <distribution.h>


/**********************************************/
//	distribution CLASS
/**********************************************/


//-------------------------------------------- constructors
distribution::distribution(double dx0) {
	// create dummy element
	key k; k.dummy=true;
	element_iterator neighbour[2*DIMS];
	elements.insert(dist_pair(k, element(elements.end()))); // new element for dummy
	element_iterator it=elements.begin();
	for (int n=0;n<DIMS;n++) {
		it->second.neighbour[LEFT] =elements.begin();	//	to left
		it->second.neighbour[RIGHT]=elements.begin();	//	to right
	}
	it->second.p=0.0;
	dx=dx0;
}

distribution::distribution() {
	distribution(0);
}


//-------------------------------------------- at -- O(log n)
element_iterator distribution::at(key k) {
	element_iterator it = elements.find(k);
	return it==elements.end() ? elements.begin() : it;	// else return dummy element
}

//-------------------------------------------- initialise distribution with a Gaussian
void distribution::initialise(double x[DIMS], double variance) {	//	Gaussian point initialise
	key k;
	for (int n=0;n<DIMS;n++) {
		k[n]=(int)(x[n]/dx - 0.5);
	}
	// single point
	insert(k,1);
	// then blur it
	int N=(int)(9*variance/(4*dx*dx));
	for (int n=0;n<N;n++) {
		update(0);
		blur(variance/(double)N);
	}
}

//-------------------------------------------- test_case
// Initialise distribution with LeVeque's test case from example 20.1
void distribution::test_case() {
	assert(DIMS==2);
	//            { 1 if 0.1 < x < 0.6 and -0.25 < y < 0.25
	// q(x,y,0) = { 1-r/0.35 if r:=sqrt((x+0.45)^2 + y^2) < 0.35
	//            { 0 otherwise
	key k;
	for (double x=-1.0; x<1.0; x+=dx) {
		for (double y=-1.0; y<1.0; y+=dx) {
			double r=sqrt( (x+0.45)*(x+0.45) + y*y);
			if ( (x>0.1) && (x<0.6) && (y>-0.25) && (y<0.25) ) {
				// the cuboid
				k[0]=(int)(x/dx - 0.5);
				k[1]=(int)(y/dx - 0.5);
				insert(k,0);
				at(k)->second.p=1; // set value, rather than add it
			} else if (r<0.35) {
				// the cone
				k[0]=(int)(x/dx - 0.5);
				k[1]=(int)(y/dx - 0.5);
				insert(k,0);
				at(k)->second.p=1-r/0.35;
			}
		}
	}
	update(0);
}

//-------------------------------------------- insert
element_iterator distribution::insert(key k, double p) {
	// if k is new, insert and update all neighbours
	
	pair<element_iterator, bool> newpair;
	// insert method leaves alone if element already exists
	newpair = elements.insert(dist_pair(k,element(elements.begin())));	//	O(log n)
	element_iterator it = newpair.first;

	if (newpair.second) {	// if element is new
		//	update element with neighbourhood
		for (int n=0;n<DIMS;n++) {
			// create a key to appropriate neighbour
			// to left
			it->second.neighbour[LEFT] =update_neighbour(k.left(n), it, RIGHT);	//	update neighbour
			//	to right
			it->second.neighbour[RIGHT]=update_neighbour(k.right(n), it, LEFT);	//	update neighbour
		}
		// the new element's vels
		it->second.p = p;
		for (int n=0;n<DIMS;n++) {									//	insert vel at right edge of box
			double edge[DIMS];
			for (int m=0;m<DIMS;m++) edge[m]=x(k[m]) + (m==n)*dx/2;	//	position of box edge
			it->second.v[n]=u(edge,n);
			it->second.w[n]=it->second.v[n]*(it->second.v[n] > 0.0);
			it->second.u[n]=it->second.v[n]*(it->second.v[n] < 0.0);
		}
	}
	else {	// element exists, so add p
		it->second.p += p;
	}
	return it;
}

//-------------------------------------------- update_neighbour -- O(log n)
inline element_iterator distribution::update_neighbour(const key &k, element_iterator caller, int d) {
	element_iterator it = elements.find(k);
	it==elements.end() ? it=elements.begin() : it->second.neighbour[d]=caller;
	return it;
}

//-------------------------------------------- save distribution to a text file
void distribution::save(string fname, double t) {
	element_iterator it;
	key k;
	stringstream os;
	os << fixed << setprecision(2) << t << "_" << fname; // convert time to string
	fname = os.str();
	ofstream file;
	file.open(fname.c_str(), ios::trunc);
	it=elements.begin(); ++it;	//	skip dummy element
	while (it != elements.end()) {
		for (int n=0;n<DIMS;n++) {
			file << x(it->first.ijk[n]) << "\t";	//	x
		}
		for (int n=0;n<DIMS;n++) {
			file << it->first.ijk[n] << "\t";		//	ijk
		}
		for (int n=0;n<DIMS;n++) {
			file << it->second.v[n] << "\t";		//	v
		}
		file << it->second.p << endl;				//	p
		++it;
	}
	file.close();
}


//-------------------------------------------- 
//	puts blur of p into p -- 1 step
inline bool distribution::blur(double variance) {
	//	Gaussian blur
	element_iterator it;
	double c = exp(-dx*dx/(2.0*variance));
	G_here = 1.0/(1.0+2.0*c);
	G_neighbour = (1.0-G_here)/2.0;
	//	repeat 1D in each dimension, rather than do Nd kernel (using superposition)
	//	just look at neighbour -- assume second-nearest contribution is small
	for (int n=0;n<DIMS;n++) {
		for (it=elements.begin()++; it != elements.end(); ++it) {
			it->second.temp = it->second.p;
		}
		for (it=elements.begin()++; it != elements.end(); ++it) {
			it->second.p = G_here*it->second.temp;
			it->second.p += G_neighbour*(it->second.neighbour[LEFT]->second.temp);
			it->second.p += G_neighbour*(it->second.neighbour[RIGHT]->second.temp);
		}
	}
	elements.begin()->second.p = 0.0;
	return true;
}

//-------------------------------------------- update -- O(n log n) -- the bottleneck
void distribution::update(double eps) {
	// look through list
	element_iterator it=elements.begin()++;
	while (it != elements.end()) {	//	O(n)
		for (int n=0;n<DIMS;n++) {
			// make sure neighbours exist
			//	update element with neighbourhood
			// if element is new,
			if (it->second.neighbour[LEFT]->first.dummy && it->second.p > eps) {
				//	look to left
				insert(it->first.left(n), 0);
			}
			if (it->second.neighbour[RIGHT]->first.dummy && it->second.p > eps) {
				//	look to left
				insert(it->first.right(n), 0);
			}
		}
		if (it->second.p < eps && !(big_neighbours(it, eps))) {
			element_iterator temp=it;
			temp++;
			drop(it);
			it=temp;
		} else ++it;
		//	if no big neighbours, drop
		//	otherwise it's on periphery
	}
	// normalise
	//*this *= 1.0/sum();
}

//-------------------------------------------- drop (by key) -- O(log n)
int distribution::drop(key k) {
	return drop(at(k));
}

//-------------------------------------------- drop (by iterator) -- O(1)
inline int distribution::drop(element_iterator it) {
	if (it->first.dummy) return 1;	//	never drop the dummy element
	// it points to the element to be dropped
	// update neighbourhood to point to dummy
	element_iterator dummy=elements.begin();
	for (int n=0;n<DIMS;n++) {
		// neighbours must point to dummy before dropping
		// to left
		//	update neighbour
		it->second.neighbour[RIGHT]->second.neighbour[LEFT] = dummy;
		it->second.neighbour[LEFT]->second.neighbour[RIGHT] = dummy;
	}
	elements.erase(it);	//	drop element
	return 0;
}

//-------------------------------------------- big_neighbours
inline bool distribution::big_neighbours(element_iterator it, double eps) {
	// loop over neighbours, if it or any has p > eps, return true
	if (it->second.p > eps) return true;	//	is itself big
	for (int n=0;n<DIMS;n++) {
		if (((it->second.neighbour[LEFT])->second.p) > eps) {
			//			cout << " has a big neighbour " << endl;
			return true;
		}
		if (((it->second.neighbour[RIGHT])->second.p) > eps) {
			//			cout << " has a big neighbour " << endl;
			return true;
		}
	}
	//	cout << " no big neighbour " << endl;
	return (it->first.dummy);
}

//-------------------------------------------- assignment
distribution& distribution::operator=(const distribution& Q) {
	if (this != &Q) {		// can't assign to self
		dx=Q.dx;
		elements=Q.elements;
	}
	dx=Q.dx;
	return *this;			//	return reference for multiple assignments
}

//-------------------------------------------- addition
distribution distribution::operator+(distribution& Q) {
	distribution S;
	S=*this; // S inherits everything
	S+=Q;
	return S;
}

//-------------------------------------------- addition +=
distribution& distribution::operator+=(distribution& Q) {
	element_iterator it;
	// loop over Q
	for (it=(Q.elements).begin() ; it != (Q.elements).end(); ++it ){
		insert(it->first, it->second.p + (elements[it->first]).p);	//	Pijk=Pijk + Qijk
	}
	return *this;
}

//-------------------------------------------- subtraction
distribution distribution::operator-(distribution& Q) {
	distribution S;
	S=*this; // S inherits everything
	S-=Q;
	return S;
}

//-------------------------------------------- subtraction -=
distribution& distribution::operator-=(distribution& Q) {
	element_iterator it;
	// loop over Q
	for (it=(Q.elements).begin() ; it != (Q.elements).end(); ++it ){
		insert(it->first, -it->second.p + (elements[it->first]).p);	//	Pijk=Pijk - Qijk
		
	}
	return *this;
}

//-------------------------------------------- multiplication
distribution distribution::operator*(double c) {
	distribution S;
	S=*this; // S inherits everything
	S*=c;
	return S;
}

//-------------------------------------------- multiplication *=
distribution& distribution::operator*=(double c) {
	element_iterator it;
	for (it=elements.begin()++ ; it != elements.end(); ++it ){
		it->second.p *= c;
	}
	return *this;
}

//-------------------------------------------- sum over elements
double distribution::sum() {
	element_iterator it;
	double psum=0;
	for (it=elements.begin()++ ; it != elements.end(); ++it ){
		psum+=it->second.p;
	}
	return psum;
}

//-------------------------------------------- max v over elements
double distribution::max_v() {
	element_iterator it;
	double maxv=0;
	for (it=elements.begin()++ ; it != elements.end(); ++it ){
		double mv=0;
		for (int n=0;n<DIMS;n++) {
			mv+=(it->second.v[n])*(it->second.v[n]);
		}
		maxv=max(sqrt(mv),maxv);
	}
	return maxv;
}

//-------------------------------------------- return size
int distribution::size() {
	return elements.size();
}

//--------------------------------------------
// Kullback-Liebler: D_KL(P||K)
// returns expected bits if message coded using Q where P is true distribution
// http://en.wikipedia.org/wiki/Kullback–Leibler_divergence
double distribution::D_KL(distribution Q) {
	// P and Q required to sum to 1
	double p_sum = this->sum();
	double q_sum = Q.sum();
	double H=0.0;
	element_iterator it, it2;
	assert(dx==Q.dx); // same grid, else interpolation would be necessary
	for (it=elements.begin()++ ; it != elements.end(); ++it ){
		double p = it->second.p/p_sum;
		it2 = Q.at(it->first); // look up same position in Q as in P, require Q defined everywhere p>0
		double q = it2->second.p/q_sum;
		H += (p > 0.0) ? p*log(p/q) : 0.0;	// 0*log(0) is defined as 0
	}
	return H/log(2);
}

//--------------------------------------------
// Kullback-Liebler divergence relative to static PDF m
// See: http://en.wikipedia.org/wiki/Kullback–Leibler_divergence
// H is information (in nats) going from PDF m to PDF p

// See http://en.wikipedia.org/wiki/Limiting_density_of_discrete_points

// from http://en.wikipedia.org/wiki/Maximum_entropy_thermodynamics :
/* The simple definition of Shannon entropy ceases to be directly applicable for
 random variables with continuous probability distribution functions. Instead
 the appropriate quantity to maximise is the "relative information entropy," 
 Hc is the negative of the Kullback-Leibler divergence, or discrimination
 information, of m(x) from p(x), where m(x) is a prior invariant measure for the
 variable(s). The relative entropy Hc is always less than zero, and can be
 thought of as (the negative of) the number of bits of uncertainty lost by
 fixing on p(x) rather than m(x). Unlike the Shannon entropy, the relative
 entropy Hc has the advantage of remaining finite and well-defined for
 continuous x, and invariant under 1-to-1 coordinate transformations. The two
 expressions coincide for discrete probability distributions, if one can make
 the assumption that m(xi) is uniform - i.e. the principle of equal a-priori
 probability, which underlies statistical thermodynamics.*/

// Kullback-Liebler of P vs Gaussian
double distribution::H(double var) {
	*this *= 1.0/sum();
	double H=0.0;
	element_iterator it;
	double m, c;
	double dV = pow(dx,DIMS);
	// c is normalising factor for multivariate Gaussian
	c=pow(2.0*M_PI,(double)DIMS/2.0);
	for (it=elements.begin()++ ; it != elements.end(); ++it ){
		// multivariate Gaussian
		double dsq=0.0;
		for (int n=0; n<DIMS; n++) dsq+=x(it->first.ijk[n])*x(it->first.ijk[n]);
		m = dV*exp(-dsq/(2*var))/(sqrt(var)*c);
		H += (it->second.p > 0.0) ? it->second.p*log(it->second.p/m) : 0.0;
	}
	return H;
}

// Shannon-Jaynes entropy
double distribution::H() {
	*this *= 1.0/sum();
	double H=0.0;
	element_iterator it;
	double m;
	double dV = pow(dx,DIMS);
	m=dV*1;
	for (it=elements.begin()++ ; it != elements.end(); ++it ){
		H -= (it->second.p > 0.0) ? it->second.p*log(it->second.p/m) : 0.0;
	}
	return H;
}

//-------------------------------------------- 
// Kullback-Liebler divergence relative to p at previous timestep
// Doesn't really make sense since timestep dependent because
// H does not obey triangle inequality
double distribution::Hdot(double dt) {
	double Hdot=0.0;
	double dV = pow(dx,DIMS);
	element_iterator it;
	for (it=elements.begin()++ ; it != elements.end(); ++it ){
		Hdot += (it->second.p > 0.0) ? it->second.p*log(it->second.p/it->second.temp) : 0.0;
	}
	return Hdot/dt;
}

//-------------------------------------------- RHS or dP/dt
int distribution::step(double dt) {
	element_iterator it;
	element_iterator i;
	element_iterator j;
	element_iterator k;
	element_iterator m;
	
	//	DRIFT FLUX / ADVECTION TERMS
	for (it=elements.begin()++; it != elements.end(); ++it) {
		for (int n=0;n<DIMS;n++) {
			it->second.flux[n] = (it->second.w[n])*(it->second.p) + (it->second.u[n])*(it->second.neighbour[RIGHT]->second.p);
		}
	}

	for (it=elements.begin()++; it != elements.end(); ++it) {
		for (int n=0;n<DIMS;n++) {
			i=it->second.neighbour[LEFT];
			double f = (it->second.p - i->second.p)*dt/(2*dx);
			
			// Corner transport upwind flux terms
			for (int e=0;e<DIMS;e++) {
				bool comp;
				comp = (e!=n);
				it->second.flux[e] -= comp*(it->second.w[e] * i->second.w[n] * f);
				j = it->second.neighbour[LEFT];
				j->second.flux[e] -= comp*(j->second.u[e] * i->second.w[n] * f);
				i->second.flux[e] -= comp*(i->second.u[e] * i->second.w[n] * f);
				m = i->second.neighbour[LEFT];
				m->second.flux[e] -= comp*(m->second.u[e] * i->second.u[n] * f);
			}

			// Second-order correction flux term
			double theta, t;
			theta = (i->second.v[n] > 0) ?	
				(i->second.p - i->second.neighbour[LEFT]->second.p)/(it->second.p - i->second.p) :
				(it->second.neighbour[RIGHT]->second.p - it->second.p)/(it->second.p - i->second.p);
			t = abs(i->second.v[n]);
			i->second.flux[n] += t*(dx/dt - t)*f*mc(theta);
		}
	}

	//	DIFFUSION TERMS
	VAR ? blur(VAR*dt) : 0;	//	just blurs, so p -> p + diffusion term
	
	//	BOTH TERMS ADDED, EXPLICIT EULER INTEGRATION IN TIME
	//	so p + diffusion -> p + diffn + advection
	for (it=elements.begin()++; it != elements.end(); ++it) {
		for (int n=0;n<DIMS;n++) {
			//	march, explicit Euler for advection term
			it->second.p -= (it->second.flux[n] - it->second.neighbour[LEFT]->second.flux[n])*dt/dx;
		}
	}
	
	//*this *= 1.0/sum();
	return 0;
}

//-------------------------------------------- flux limiter functions
inline double distribution::mc(double theta) {
	// monotonised central difference limiter
	return max(0.0, min((1.0+theta)/2.0, min(2.0, 2.0*theta)));
}

inline double distribution::vl(double theta) {
	// van Leer limiter
	return min(0.0, (theta + abs(theta))/(1.0+abs(theta)));
}


//-------------------------------------------- Bayesian update
//	observation is obs, with variance var, on dimension n only
int distribution::observe(double obs, double var, int n) {
	var *= 2.0;
	for (element_iterator it=elements.begin()++; it != elements.end(); ++it) {
		double d = (x(it->first.ijk[n])-obs);
		it->second.p *= exp(-d*d/var);
	}
	*this *= 1.0/sum();
	return 0;
}
