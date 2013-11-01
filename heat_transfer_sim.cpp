
#include <iostream>
#include <complex>
#include <cstdlib>
#include "fieldFunctions.hpp"
#include "SMatrix.hpp"
#include "mlgeo.hpp"
#include "materials.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);

double fluxKpInt(const mlgeo &g, const int *l, const double *zl, const double *nHat,
				 const int *s, int Nl, int Ns, double k0);
double fluxKpWInt(const epsfn *eps, const double *d, const int *l, const double *zl,
				  const double *nHat, const int *s, int numLayers, int Nl, int Ns, 
				  double lscale, double temp);

const double c0 = 299798452;

int main() {
	// frequency and wavevector
	double w1 = 1.0e14; 
	double w2 = 1.0e15;
	int Nw = 200;

	double k1 = 0.01;
	double k2 = 45;
	int Nk = 150;
	
	const bool omegaK = false;
	const bool omega = true;
	const bool print_eps = true;
	//const bool boltzmann = true;
	const double temp = 1000.;

	const double lscale = 1e-6; // nm
	
	// pararmeters for two multilayer stacks 
	// of alternating materials/thicknesses, with some thickness
	// separating them.  Vacuum on either side
	//   eps0 ||  ... || eps2,d2 || eps1,d1 || epsInt,di || eps1,d1 || eps2,d2 || ... || epsN
	const int numLayersPerMat = 2;
	const epsfn eps0 = epsVac;
	const epsfn eps1 = epsCr;
	const epsfn eps2 = [=] (double w)  { return epsConst(w,3.9); }; // SiO2=3.9
	const epsfn epsN = epsVac;
	const epsfn epsInt = epsVac; // spacer layer
	const double ff = 0.2;
	const double totalD = 0.1; 
	const double d1 = ff * totalD; // thickness of eps1 (um)
	const double d2 = (1.-ff) * totalD;
	const double di = 0.5; // distance between stacks
	const int numLayers = numLayersPerMat*2 + 3;
	const int intLayer = numLayersPerMat + 1;
	const int lastLayer = 2 * intLayer;

	epsfn *eps = new epsfn[numLayers];
	cdouble *epsV = new cdouble[numLayers];
	double *d = new double[numLayers-2];
	const int numEmitters = numLayersPerMat;

	// arrays with eps functions & layer thicknesses
	eps[0] = eps0;
	eps[intLayer] = epsInt;
	eps[lastLayer] = epsN;
	d[intLayer-1] = di;
	for(int i=0; i<numLayersPerMat; ++i) {
		eps[intLayer-1-i] = (i%2==0) ? eps1 : eps2;
		eps[intLayer+1+i] = (i%2==0) ? eps1 : eps2;
		d[intLayer-2-i] = (i%2==0) ? d1 : d2;
		d[intLayer+i] = (i%2==0) ? d1 : d2;
	}

	if(print_eps) {
		for(int i=0; i<numLayers; ++i)
			epsV[i] = eps[i](w1);
		mlgeo g = mlgeo(epsV, d, numLayers-1);
		g.print();
		printEps(eps1, w1, w2, Nw);
		return 0;
	}

	// compute emitting layers
	int *s = new int[numEmitters];
	for (int i = 0; i < numEmitters; ++i) 
		s[i] = i + 1;

	if(omegaK || omega) {
		for (int iw = 0; iw < Nw; ++iw) {
			double w = (Nw != 1) ?  w1 + iw * (w2 - w1) / (Nw - 1) : w1;
			double k0 = w * lscale / c0;
						
			for (int i=0; i<numLayers; ++i)
				epsV[i] = eps[i](w);
			mlgeo g = mlgeo(epsV, d, numLayers-1);
		
			if (omegaK) {
				for (int ik = 0; ik < Nk; ++ik) {
					double kp = (Nk!=1) ? (k1 + ik*(k2-k1)/(Nk-1))*k0 : k1*k0;
					SMatrix S = SMatrix(g, k0, kp);
					double q1 = flux(g, S, intLayer, di/2., s, numLayersPerMat, 1.);
					double q2 = flux(g, S, lastLayer, di/2., s, numLayersPerMat, 1.);
					std::cout << w << " " << kp/k0 << " " << (q1 - q2) << std::endl;
				}
			} else {
				int l[2] = {intLayer, lastLayer};
				double zl[2] = {di/2., di/2.};
				double nHat[2] = {1.,-1.};
				double q = fluxKpInt(g, l, zl, nHat, s, 2, numEmitters, k0);
				std::cout << w << " " << q << std::endl; 
				//std::cout << w << " " << q / flux_bb(k0) << std::endl; 
			}
		}
	} else {
		int l[2] = {intLayer, lastLayer};
		double zl[2] = {di/2., di/2.};
		double nHat[2] = {1.,-1.};
		double q = fluxKpWInt(eps, d, l, zl, nHat, s, numLayers, 2, numEmitters, lscale, temp);
		std::cout << di << " " << q / flux_bb_int(lscale, temp) << std::endl;
	}
	delete[] s;
	delete[] d;
	delete[] eps;
	delete[] epsV;
}
