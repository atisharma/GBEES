/*		A S Sharma June 2009		*/

#ifndef ELEMENT_H
#define ELEMENT_H

#include <defs.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>

using namespace std;

#define ADD 1
#define LEFT 2*n
#define RIGHT 2*n+1

/**********************************************/
//	key CLASS
/**********************************************/
class key {
private:
public:
	key();
	key(int argkey[]);
	int ijk [DIMS];
	bool dummy;
	inline key left(int) const;
	inline key right(int) const;
	inline int& operator[](const int);
	inline key& operator=(const key&);
	inline key& operator=(const int argkey[]);
	inline friend bool operator<(const key&, const key&);
	inline friend bool operator>(const key&, const key&);
	inline friend bool operator==(const key&, const key&);
	inline friend bool operator!=(const key&, const key&);
	void display();
	void display() const;
};


//--------------------------------------------
inline key key::left(int i) const {	
	key k;
	k.dummy=false;
	for (int n=0;n<DIMS;n++) {
		k.ijk[n]=ijk[n];
	}
	k.ijk[i]=ijk[i]-1;
	return k;
}

//--------------------------------------------
inline key key::right(int i) const {	
	key k;
	k.dummy=false;
	for (int n=0;n<DIMS;n++) {
		k.ijk[n]=ijk[n];
	}
	k.ijk[i]=ijk[i]+1;
	return k;
}

//--------------------------------------------
inline int& key::operator[](const int i) {
	//assert(0 <= i && i < DIMS);
	return ijk[i];
}

//--------------------------------------------
inline key& key::operator=(const key& v) {
	for (int n=0;n<DIMS;n++) {
		ijk[n]=v.ijk[n];
	}
	return *this;
}

//--------------------------------------------
inline key& key::operator=(const int argkey[]) {
	for (int n=0;n<DIMS;n++) {
		ijk[n]=argkey[n];
	}
	return *this;
}

//-------------------------------------------- for sorting
inline bool operator<(const key& a, const key& b) {
	int c=2;
	for (int n=0;n<DIMS;n++) {
		(a.ijk[n]<b.ijk[n]) ? c=1 : NULL;
		(a.ijk[n]>b.ijk[n]) ? c=0 : NULL;
	}
	(a.dummy > b.dummy) ? c=1 : NULL;
	(a.dummy < b.dummy) ? c=0 : NULL;
	(c==2) ? c=0 : NULL;
	return c;
}

//--------------------------------------------
inline bool operator>(const key& a, const key& b) {
	return ((a==b)||(a<b)) ? 0 : 1;
}

//--------------------------------------------
inline bool operator==(const key& a, const key& b) {
	return !(a<b) && !(b<a);
}

//--------------------------------------------
inline bool operator!=(const key &a, const key& b) {
	return !(a==b);
}


/**********************************************/
//	element CLASS
/**********************************************/
class element {
private:
	typedef map<const key, element>::iterator element_iterator;
public:
	double p;
	void display();
	//element(element_iterator[2*DIMS], double, double[DIMS]);
	element();
	element(element_iterator);
	element_iterator neighbour[2*DIMS];	//	pointers to neighbours
	double v[DIMS];						//	velocities on boundary
	double u[DIMS];						//	velocities on boundary
	double w[DIMS];						//	velocities on boundary
	double flux[DIMS];					//	fluxes through boundary
	double temp;						//	used in blurring; a temp variable
};

#endif
