
#include <iostream>
#include <complex>
#include "fieldFunctions.hpp"
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
	const int numLayersPerMat = 4;
	const epsfn eps0 = epsVac;
    //const epsfn eps1 = [=] (double w) { return epsSiDoped(w,5.e19); };
	const epsfn eps1 = epsSiC;
	const epsfn eps2 = [=] (double w)  { return epsConst(w,3.9); }; // SiO2=3.9
	const epsfn epsN = epsVac;
	const epsfn epsInt = epsVac; // spacer layer
	const double ff = 0.2;
	const double totalD = 250.;
	const double d1 = ff * totalD; // thickness of eps1 (um)
	const double d2 = (1.-ff) * totalD;
	const double di = 100; // distance between stacks
	const double units = 1e-9; // nm

	bool boltzmann = false;
	bool omegaK = true;
	bool print_eps = false;

	// frequency and wavevector
	double w1 = 1.0e14; 
	double w2 = 2.4e14;
	int Nw = 100;

	double k1 = 1.01;
	double k2 = 40;
	int Nk = 100;

	int numLayers = numLayersPerMat*2 + 3;
	epsfn *eps = new epsfn[numLayers];
	cdouble *epsV = new cdouble[numLayers];
	double *d = new double[numLayers-2];

	// setup stack, eps & d values
	int intLayer = numLayersPerMat + 1;
	int lastLayer = 2 * intLayer;
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
		mlgeo *g = new mlgeo(epsV, d, numLayers-1);
		g->print();
		printEps(eps1, w1, w2, Nw);
		return 0;
	}

	double w, k0, q1, q2, kp, norm = 1.;
	for (int iw=0; iw<Nw; ++iw) {
		w = (Nw!=1) ?  w1 + iw*(w2-w1)/(Nw-1) : w1;
		k0 = w/3.e8 * units;
		
		if(boltzmann)
			norm = meanEnergy(w,300);
					
		for(int i=0; i<numLayers; ++i)
			epsV[i] = eps[i](w);
		mlgeo *g = new mlgeo(epsV, d, numLayers-1);
	
		if(omegaK) {
			for(int ik=0; ik<Nk; ++ik) {
				kp = (Nk!=1) ? (k1 + ik*(k2-k1)/(Nk-1))*k0 : k1*k0;
				q1 = 0; q2 = 0;
				SMatrix S = SMatrix(g, k0, kp);
				/*for(int s = 1; s <= numLayersPerMat; ++s) {
					q1 += flux(g, intLayer, s, di/2., k0, kp, 1.);
					q2 += flux(g, lastLayer, s, di/2., k0, kp, 1.);
				}
				std::cout << w << " " << kp/k0 << " " << norm*(q1-q2)/units << std::endl;
				*/
				int s = 1;
				double zs = g->d(s) / 2.;
				pwaves pTE, pTM;
				pWavesL(&S, intLayer, s, TE, pTE); 
				pWavesL(&S, intLayer, s, TM, pTM); 
				cdouble qa1 = gfFluxTE(&S, pTE, intLayer, s, di/2., zs, false)
							+ gfFluxTM(&S, pTM, intLayer, s, di/2., zs, false);
				cdouble qa2 = gfFluxSP(&S, pTE, pTM, intLayer, s, di/2., zs);
				std::cout << w << " " << kp/k0 << " " << qa1 << " " << qa2 << std::endl;
			}
		} else {
			q1 = 0; q2 = 0;
			for(int s = 1; s <= numLayersPerMat; ++s) {
				q1 += fluxKpInt(g, intLayer, s, di/2., k0, 1.);
				q2 += fluxKpInt(g, lastLayer, s, di/2., k0, 1.);
			}
			std::cout << w << " " << norm*(q1-q2)/(units*units) << std::endl; // e.g. \um^2 -> \m^2
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
