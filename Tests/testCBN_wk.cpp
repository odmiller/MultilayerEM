
#include <iostream>
#include <complex>
#include "multilayer.hpp"
#include "SMatrix.hpp"
#include "materials.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);

int main() {
	// pararmeters for two multilayer stacks 
	// of alternating materials/thicknesses, with some thickness
	// separating them.  Vacuum on either side
	//   eps0 || eps1,d1 || eps2,d2 || ... || epsInt,di || ... || eps1,d1 || epsN
	const int numLayersPerMat = 1;  // odd, for mirror symmetry
	const epsfn eps0 = epsVac;
    const epsfn eps1 = epsCBN;
	const epsfn eps2 = epsVac; // [=] (double d)  { return epsConst(d,3.9); }; // SiO2=3.9
	const epsfn epsN = epsVac;
	const epsfn epsInt = epsVac; // spacer layer
	const double d1 = 2; // thickness of eps1 (um)
	const double d2 = 0;
	const double di = 0.1; // distance between stacks

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
	double w, k0, kp, w1, w2, Nw, k1, k2, Nk, q, t;
	int iw, ik;
	const double hbar = 6.62606957e-34/2/M_PI;
	const double kB = 1.3806488e-23;
	w1 = 1.6e14; 
	w2 = 2.6e14;
	k1 = 1.01;
	k2 = 30;
	Nw = 50;
	Nk = 50;

	for (iw=0; iw<Nw; ++iw) {
		w = w1 + iw*(w2-w1)/(Nw-1);
		k0 = w/3.e14; // 1/um
		
		for(int i=0; i<numLayers; ++i)
			epsV[i] = eps[i](w);
		//printEps(epsCBN, 0.8*w,1.2*w,50);

		mlgeo *g = new mlgeo;
		initGeo(epsV, d, numLayers-1, g);

		t = hbar*w / (exp(hbar*w/(kB*300)) - 1); // Boltzmann Theta

		for (ik=0; ik<Nk; ++ik) {
			kp = (k1 + ik*(k2-k1)/(Nk-1)) * k0;
			q = flux(g, 2, 1, 0.05, k0, kp, 1.);
			std::cout << w << " " << kp << " " << q*t*1e6 << std::endl; // \um -> \m
		}
		delete g;
	}

	//std::cout << "\nflux from all of s: " << h << std::endl;
	//std::cout << "t: " << t << std::endl;
	//std::cout << "flux * Boltzmann: " << h*t << std::endl;
}
