/*		A S Sharma June 2009		*/

#include <element.h>

/**********************************************/
//	key CLASS
/**********************************************/

//--------------------------------------------
key::key(){
	for (int n=0;n<DIMS;n++) { ijk[n]=0; }
	dummy=false;
}

//--------------------------------------------
key::key(int argkey[]){
	for (int n=0;n<DIMS;n++) {
		ijk[n]=argkey[n];
	}
	dummy=false;
}

//--------------------------------------------
void key::display(){
	if (!dummy) {
		cout << "(";
		for (int n=0;n<DIMS;n++) cout << "\t" << ijk[n];
		cout << ")";
		cout << endl;
		return;
	}
	cout << ".";
	return;
}	

//--------------------------------------------
void key::display() const {
	display();
}	

/**********************************************/
//	element CLASS
/**********************************************/

/*
//-------------------------------------------- constructor
element::element (element_iterator newneighbour[2*DIMS], double newp, double vel[DIMS]) {
	p=newp;
	for (int n=0;n<DIMS;n++) {
		neighbour[LEFT]=newneighbour[LEFT];		//	box to left
		neighbour[RIGHT]=newneighbour[RIGHT];	//	box to right
		v[n]=vel[n];
		w[n]=vel[n]*(vel[n] > 0.0);
		u[n]=vel[n]*(vel[n] <= 0.0);
		flux[n]=0.0;
	}
	temp=0.0;
}
*/
//-------------------------------------------- constructor
element::element (element_iterator i) {
	p=0.0;
	for (int n=0;n<DIMS;n++) {
		neighbour[LEFT]=i;
		neighbour[RIGHT]=i;
		v[n]=0.0;
		w[n]=0.0;
		u[n]=0.0;
		flux[n]=0.0;
	}
	temp=0.0;
}

element::element () {
	p=0.0;
	for (int n=0;n<DIMS;n++) {
		v[n]=0.0;
		w[n]=0.0;
		u[n]=0.0;
		flux[n]=0.0;
	}
	temp=0.0;
}

