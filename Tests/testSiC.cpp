
/* This code is used to test both the w-k image plot and the 
 * q-w plot (for different separations) in the paper: 
 * Francoeur et. al. J. Phys. D. 43, 075501 (2010),
 * which is referred to as Ref2 in some of the output data
 */

#include <iostream>
#include <complex>
#include "multilayer.hpp"
#include "SMatrix.hpp"
#include "mlgeo.hpp"
#include "materials.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);

double fluxKpInt(const mlgeo *g, int l, int s, double zl, double k0, double nHat);

int main() {
	// pararmeters for two multilayer stacks 
	// of alternating materials/thicknesses, with some thickness
	// separating them.  Vacuum on either side
	//   eps0 || eps1,d1 || eps2,d2 || ... || epsInt,di || ... || eps1,d1 || epsN
	const int numLayersPerMat = 1;  // odd, for mirror symmetry
	const epsfn eps0 = epsVac;
    const epsfn eps1 = epsSiC;
	const epsfn eps2 = epsVac; // [=] (double d)  { return epsConst(d,3.9); }; // SiO2=3.9
	const epsfn epsN = epsVac;
	const epsfn epsInt = epsVac; // spacer layer
	const double d1 = 10; // thickness of eps1 (um)
	const double d2 = 0;
	const double di = 10; // distance between stacks
	const double units = 1e-9; // nm
	const int s = 1; // source layer

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
	double w, k0, w1, w2, Nw, q1, q2, t;
	int iw;
	const double hbar = 6.62606957e-34/2/M_PI;
	const double kB = 1.3806488e-23;
	w1 = 1.5e14; 
	w2 = 1.9e14;
	Nw = 100;

	double kp, k1, k2;
	int Nk, ik;
	k1 = 1.01;
	k2 = 200;
	Nk = 100;

	//printEps(epsSiC, w1, w2, Nw);
	bool omegaK = true;
	for (iw=0; iw<Nw; ++iw) {
		w = (Nw!=1) ?  w1 + iw*(w2-w1)/(Nw-1) : w1;
		k0 = w/3.e8 * units; // 1/um
		
		for(int i=0; i<numLayers; ++i)
			epsV[i] = eps[i](w);

		mlgeo *g = new mlgeo(epsV, d, numLayers-1);

		t = hbar*w / (exp(hbar*w/(kB*300)) - 1); // Boltzmann Theta
		
		if(omegaK) {
			for(ik=0; ik<Nk; ++ik) {
				kp = (Nk!=1) ? (k1 + ik*(k2-k1)/(Nk-1))*k0 : k1*k0;

				q1 = flux(g, 2, s, di/2., k0, kp, 1., 5);
				q2 = flux(g, 4, s, di/2., k0, kp, 1.);
				std::cout << w << " " << kp/k0 << " " << (q1-q2)*t/units << std::endl;
			}
		} else {
			q1 = fluxKpInt(g, 2, 1, di/2., k0, 1.);
			q2 = fluxKpInt(g, 4, 1, di/2., k0, 1.);
			std::cout << w << " " << (q1-q2)*t/(units*units) << std::endl; // e.g. \um^2 -> \m^2
		}
		delete g;
	}

	delete[] d;
	delete[] eps;
	delete[] epsV;
	//std::cout << "\nflux from all of s: " << h << std::endl;
	//std::cout << "t: " << t << std::endl;
	//std::cout << "flux * Boltzmann: " << h*t << std::endl;
}
