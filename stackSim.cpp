
#include <iostream>
#include <complex>
#include "multilayer.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);

int main() {
	//epsfn geoEps[] = { epsVac, epsAu, epsAg, epsVac };	
	
	// pararmeters for two multilayer stacks 
	// of alternating materials/thicknesses, with some thickness
	// separating them.  Vacuum on either side
	//   eps0 || eps1,t1 || eps2,t2 || ... || epsInt,ti || ... || eps1,t1 || epsN
	const int numLayersPerMat = 3;  // odd, for mirror symmetry
	const epsfn eps0 = epsVac;
    const epsfn eps1 = epsAu;
	const epsfn eps2 = [=] (double d)  { return epsConst(d,3.9); }; // SiO2=3.9
	const epsfn epsN = epsVac;
	const epsfn epsInt = epsVac; // spacer layer
	const double t1 = 0.010; // thickness of eps1
	const double t2 = 0.090;
	const double ti = 0.100; // distance between stacks

	int numLayers = numLayersPerMat*2 + 3;
	epsfn *eps = new epsfn[numLayers];
	double *t = new double[numLayers-2];

	// setup stack, eps & t values
	eps[0] = eps0;
	eps[numLayers-1] = epsN;
	for(int i=0; i<numLayersPerMat; ++i) {
		eps[i+1] = (i%2==0) ? eps1 : eps2;
		eps[numLayers-1-i-1] = eps[i+1];
		t[i] = (i%2==0) ? t1 : t2;
		t[numLayers-i-3] = t[i];
	}
	eps[numLayersPerMat+1] = epsInt;
	t[numLayersPerMat] = ti;

	for(int i=0; i<numLayers; ++i) {
		std::cout << eps[i](2*M_PI*3e8/500e-9) << std::endl;
		if(i<numLayers-2) 
			std::cout << t[i] << std::endl;
	}

}
