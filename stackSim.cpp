
#include <iostream>
#include <complex>
#include "multilayer.hpp"
#include "SMatrix.hpp"
#include "materials.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);

// forward decl to avoid fluxIntegrals header
double fluxZInt(mlgeo *g, int l, int s, double zl, double k0, double kp, double nHat);

int main() {
	//epsfn geoEps[] = { epsVac, epsAu, epsAg, epsVac };	
	
	// pararmeters for two multilayer stacks 
	// of alternating materials/thicknesses, with some thickness
	// separating them.  Vacuum on either side
	//   eps0 || eps1,d1 || eps2,d2 || ... || epsInt,di || ... || eps1,d1 || epsN
	const int numLayersPerMat = 1;  // odd, for mirror symmetry
	const epsfn eps0 = epsVac;
    const epsfn eps1 = epsCBN;
	const epsfn eps2 = epsVac; //[=] (double d)  { return epsConst(d,3.9); }; // SiO2=3.9
	const epsfn epsN = epsVac;
	const epsfn epsInt = epsVac; // spacer layer
	const double d1 = 0.010; // thickness of eps1 (um)
	const double d2 = 0.090;
	const double di = 0.100; // distance between stacks

	int numLayers = numLayersPerMat*2 + 3;
	epsfn *eps = new epsfn[numLayers];
	cdouble *epsV = new cdouble[numLayers];
	double *d = new double[numLayers-2];

	// setup stack, eps & d values
	eps[0] = eps0;
	eps[numLayers-1] = epsN;
	for(int i=0; i<numLayersPerMat; ++i) {
		eps[i+1] = (i%2==0) ? eps1 : eps2;
		eps[numLayers-1-i-1] = eps[i+1];
		d[i] = (i%2==0) ? d1 : d2;
		d[numLayers-i-3] = d[i];
	}
	eps[numLayersPerMat+1] = epsInt;
	d[numLayersPerMat] = di;

	// pick frequency and wavevector
	double w, k0, kp;
	w = 2*M_PI*3.e8/500e-9;
	k0 = w/3.e14; // 1/um
	//kp = k0 * sin(M_PI/6);
	kp = k0 * 3;
	for(int i=0; i<numLayers; ++i)
		epsV[i] = eps[i](w);

	mlgeo *g = new mlgeo;
	initGeo(epsV, d, numLayers-1, g);
	printGeo(g);

	SMatrix *S = new SMatrix(g, k0, kp);
	//S->printSMatrix();

	pwaves p0e, p2e, p0h, p2h;
	pWavesL(S, 0, 1, TE, p0e);
	pWavesL(S, 2, 1, TE, p2e);
	pWavesL(S, 0, 1, TM, p0h);
	pWavesL(S, 2, 1, TM, p2h);

/*	std::cout << "FIELDS IN LAYER 0" << std::endl;
	gfFlux(S, g, p0e, p0h, 0, 1, -0.1, 0.05);
	std::cout << "A.h0: " << p0h.Al << "  B.h0: " << p0h.Bl << "  C.h0" << p0h.Cl 
		<< "  D.h0: " << p0h.Dl << std::endl;
	std::cout << "A.e0: " << p0e.Al << "  B.e0: " << p0e.Bl << "  C.e0" << p0e.Cl 
		<< "  D.e0: " << p0e.Dl << std::endl;
	std::cout << "\n\n" << std::endl;

	std::cout << "FIELDS IN LAYER 2" << std::endl;
	gfFlux(S, g, p2e, p2h, 2, 1, 0.1, 0.05);
	std::cout << "A.h2: " << p2h.Al << "  B.h2: " << p2h.Bl << "  C.h2" << p2h.Cl 
		<< "  D.h2: " << p2h.Dl << std::endl;
	std::cout << "A.e2: " << p2e.Al << "  B.e2: " << p2e.Bl << "  C.e2" << p2e.Cl 
		<< "  D.e2: " << p2e.Dl << std::endl;
	std::cout << "\n\n" << std::endl;
*/
	double h = flux(g, 2, 3, 0.05, k0, kp, -1., 0.005);
	std::cout << "\nflux from zs=0.005: " << h << std::endl;
	
	double h2 = flux(g, 2, 1, 0.05, k0, kp, 1., 0.005);
	std::cout << "\nflux from zs=0.005: " << h2 << std::endl;

	// test numerical integration
	std::cout << "\nEntering numerical integration procedure:" << std::endl;
	fluxZInt(g, 2, 3, 0.05, k0, kp, -1.);
	
	// test analytical integral
	h2 = flux(g, 2, 3, 0.05, k0, kp, -1.);
	std::cout << "\nflux from all of s: " << h2 << std::endl;
}
