
#include <iostream>
#include <complex>
#include "materials.hpp"
#include "scatter.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);

// for given stack (eps,t,N==#interfaces), compute flux at zl (in layer l) from 
// emitter at zs (in layer s), for angular freq w and parallel wavevector kp
//   in Francoeur et. al., zl --> zc and zs -->zprime
//double flux(const epsfn *eps, const double *t, int N, int l, double zl, 
		//int s, double zs, double kp, double w) {
		
// for given stack (sTE,sTM,N==#interfaces), compute flux at zl (in layer l) from 
// emitter at zs (in layer s), for angular freq w and parallel wavevector kp
//   in Francoeur et. al., zl --> zc and zs -->zprime
double flux(const cdouble *sTE, const cdouble *sTM, int N, int l, double zl, 
		int s, double zs, double kp, double w) {
	// zl, zs are positions relative to left boundary of resp. layers
	//   zl<0 if l==0, zl>0 if l==N, zl \in [0,t[l]] if 0<l<N
	
	kzs = sqrt( epsS * k0 * k0 - kp * kp ); // ARGS?
	Splus = exp(-1.*II*kzs*zs);
	Sminus = exp(II*kzs*zs);

	// partial waves in s,0,N layers (A0=BN=C0=DN=0)
	Bs = S21(s,N,TE) * Splus /( 1 - S12(0,s,TE)*S21(s,N,TE) ); 
	As = S12(0,s,TE) * Bs;
	B0 = S22(0,s,TE) * Bs;
	AN = S11(s,N,TE) * (As + Splus);
	Cs = S12(0,s,TE) * Sminus /( 1 - S12(0,s,TE)*S21(s,N,TE) );
	Ds = S21(s,N,TE) * Cs;
	CN = S11(s,N,TE) * Cs;
	D0 = S22(0,s,TE) * (Ds + Sminus);
	if(l<s) { // flux to the left of emitter
		Bl = B0 / S22(0,l,TE);
		Al = S12(0,l,TE) * Bl;
		Dl = D0 / S22(0,l,TE);
		Cl = S12(0,l,TE) * Dl;
	} else { // flux to the right of emitter
		Al = AN / S11(l,N,TE);
		Bl = S21(l,N,TE) * Al;
		Cl = CN / S11(l,N,TE);
		Dl = S21(l,N,TE) * Cl;
	}

	return 0;
}

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
