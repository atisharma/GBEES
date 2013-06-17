/*		A S Sharma June 2009		*/

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <defs.h>
#include <element.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

typedef map<const key, element>::iterator element_iterator;
typedef pair<const key, element> dist_pair;

/**********************************************/
//	distribution CLASS
/**********************************************/
class distribution {
private:
	map<const key, element> elements;
	
	//	diffusion variables
	double G_here;
	double G_neighbour;
	
	inline element_iterator update_neighbour(const key&, element_iterator, int);
	bool big_neighbours(element_iterator, double);
	
	element_iterator insert(key, double);
	
	int drop(key);
	inline int drop(element_iterator);
	
	inline double mc(double);
	inline double vl(double);

	inline bool blur(double);
	
	
public:
	distribution();
	distribution(double);
	
	double dx;
	
	double x(int i) {return i*dx;};					//	x from i,j,k
	element_iterator at(key);						//	lookup; shouldn't be needed often
	void update(double);							//	reset element set based on eps & neighbour values
	void initialise(double[DIMS], double);
	void test_case();
	
	distribution& operator=(const distribution&);	//	assignment maintains iterators
	distribution& operator+=(distribution&);
	distribution operator+(distribution&);
	distribution& operator-=(distribution&);
	distribution operator-(distribution&);
	distribution& operator*=(double);
	distribution operator*(double);
	
	int step(double);
	
	int observe(double, double, int);
	
	double u(double[DIMS], int);
	
	void display();
	
	double H(double);	//	entropy relative to Gaussian measure (K-L divergence)
	double H();			//	entropy relative to uniform measure (Shannon-Jaynes)
	double Hdot(double);	//	entropy relative to last p distribution
	double D_KL(distribution Q); // K-L divergence from this distribution to Q
	double sum();
	double max_v();
	int size();
	
	void save(string, double);
};



#endif
