
#include <iostream>
#include <complex>
#include "fieldFunctions.hpp"
#include "SMatrix.hpp"
#include "mlgeo.hpp"
#include "materials.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);

double *dosKpInt(const mlgeo &g, const int *l, const double *zl, int Nl, double k0);

int main() {
	// frequency and wavevector
	double w1 = 1.3e14; 
	double w2 = 2.0e14;
	int Nw = 400;

	double k1 = 0.01;
	double k2 = 45;
	int Nk = 150;
	
	const double zl = -0.002;

	const bool omegaK = false;
	const bool print_eps = false;

	const double lscale = 1e-6; // um
	
	// parameters for a single multilayer stack
	const epsfn eps0 = epsVac;
	const epsfn eps1 = epsSiC;
	const epsfn eps2 = [=] (double w) { return epsConst(w,3.9); };
	const epsfn epsN = epsVac;
	const int numLayersPerMat = 1;
	const double ff = 0.2;
	const double totalD = 0.1;
	const double d1 = ff * totalD; // thickness of eps1
	const double d2 = (1.-ff) * totalD;
	const int numLayers = numLayersPerMat + 2;

	epsfn *eps = new epsfn[numLayers];
	cdouble *epsV = new cdouble[numLayers];
	double *d = new double[numLayers-2];

	// setup stack, eps & d values
	// single stack
	eps[0] = eps0;
	eps[numLayers-1] = epsN;
	for (int i = 1; i <= numLayersPerMat; ++i) {
		eps[i] = (i%2==1) ? eps1 : eps2;
		d[i-1] = (i%2==1) ? d1 : d2;
	}

	// alternative, simple specification for small stack	
	/*const int numLayers = 4;
	const epsfn eps2 = [=] (double w) { return epsConst(w,3.9); };
	epsfn eps[numLayers] = {epsVac, epsCr, epsAu, epsVac};
	cdouble epsV[numLayers];
	double d[numLayers-2] = {0.005, 0.005};
	*/

	if(print_eps) {
		for(int i=0; i<numLayers; ++i)
			epsV[i] = eps[i](w1);
		mlgeo g = mlgeo(epsV, d, numLayers-1);
		g.print();
		//printEps(eps1, w1, w2, Nw);
		return 0;
	}

	for (int iw = 0; iw < Nw; ++iw) {
		double w = (Nw != 1) ?  w1 + iw * (w2 - w1) / (Nw - 1) : w1;
		double k0 = w / 299798452. * lscale;
		double norm = 1.;
					
		for (int i=0; i<numLayers; ++i)
			epsV[i] = eps[i](w);
		mlgeo g = mlgeo(epsV, d, numLayers-1);
	
		if (omegaK) {
			for (int ik = 0; ik < Nk; ++ik) {
				double kp = (Nk!=1) ? (k1 + ik*(k2-k1)/(Nk-1))*k0 : k1*k0;
				SMatrix S = SMatrix(g, k0, kp);
				double q = dos(g, k0, kp, 0, zl);
				std::cout << w << " " << kp / k0 << " " << q * norm << std::endl;
			}
		} else {
			int l = 0;
			double *q = dosKpInt(g, &l, &zl, 1, k0);
			std::cout << w << " " << (*q) * norm << std::endl; 
		}
	}

	delete[] d;
	delete[] eps;
	delete[] epsV;
}
