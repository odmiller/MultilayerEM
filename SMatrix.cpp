
#include <iostream>
#include <complex>
#include "SMatrix.hpp"
#include "materials.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);
const cdouble II(0.0,1.0);

// const double c0 = 299792458;
// TODO: define an ``enum'' for TE/TM types...

// Fresnel coefficients (Sipe's notation)
// NOTE: kp should be normalized by k0! pass kp/k0
void reflTrans(cdouble eps0, cdouble eps1, double kp, int TM,
	cdouble &r, cdouble &t) {
	cdouble w0, w1;
	w0 = sqrt( eps0 - kp*kp );
	w1 = sqrt( eps1 - kp*kp );
	if(TM==1) {
		r = (w0*eps1 - w1*eps0) / (w0*eps1 + w1*eps0);
		t = 2.*sqrt(eps0)*sqrt(eps1)*w0 / (w0*eps1 + w1*eps0);
	} else {
		r = (w0-w1) / (w0+w1);
		t = 2.*w0 / (w0+w1);
	} 
}

// fwd decl's of some utility fxn's
void eyeS(cdouble *S, int i);
int ind11(int a, int b);
int ind12(int a, int b);
int ind21(int a, int b);
int ind22(int a, int b);

// useful for testing
void printLayerInfo(int l, cdouble eps, int pol, cdouble r, cdouble d,
					cdouble exp1, cdouble exp2, double kp, cdouble kz) {
	std::cout << "permittivity of layer " << l << " = " << eps << std::endl;
	std::cout <<"pol: " << pol << "  l: " << l << "  r: " << r << "  d: " 
		<< d << "  exp1: " << exp1 << "  exp2: " << exp2 
		<< "  kp: " << kp << "  kz: " << kz << std::endl;
}

// units of k0, kp, only matter relative to d (only k*d products ever used)
SMatrix::SMatrix(cdouble *eps, double *d, int N, double k0, double kp)
	: sTE(new cdouble[4*(2*N+1)]), sTM(new cdouble[4*(2*N+1)]) // member init. list
{
	double dl;
	cdouble kz, exp1, exp2, r, t, *s;
	int k, l; // going to set S(0,l) and S(k,N), k=0:N,l=0:N

	for(int pol=0; pol<2; ++pol) {
		s = (pol==0) ? sTE : sTM;

		// S(0,l) for l=0:N
		k = 0;
		l = 0;
		eyeS(s,ind11(k,l));
		for(l=1; l<=N; ++l) { 
			dl = (l>1) ? d[l-2] : 0;
			kz = sqrt( eps[l-1] * k0 * k0 - kp * kp );
			exp1 = exp(II*kz*dl);
			exp2 = exp(2.*II*kz*dl);
			reflTrans( eps[l-1], eps[l], kp/k0, pol, r, t );
			
			//printLayerInfo(l, eps[l], pol, r, t, exp1, exp2, kp, kz);

			s[ind11(k,l)] = ( s[ind11(k,l-1)] * t * exp1 )
					/ ( 1. - s[ind12(k,l-1)] * r * exp2 ); 
			s[ind12(k,l)] = ( s[ind12(k,l-1)] * exp2 - r )
					/ ( 1. - s[ind12(k,l-1)] * r * exp2 );
			s[ind21(k,l)] = s[ind11(k,l)] * s[ind22(k,l-1)] * r * exp1 / t + s[ind21(k,l-1)];
			s[ind22(k,l)] = s[ind22(k,l-1)] * ( r * s[ind12(k,l)] + 1. ) * exp1 / t;
		}
		
		// S(k,N) for k=N:0 (backwards recurrence relations)
		k = N;
		l = N;
		eyeS(s,ind11(k,l));
		for(k=N; k>1; --k) {  // k-1 is the new layer, from k
			dl = (k>1) ? d[k-2] : 0;
			kz = sqrt( eps[k-1] * k0 * k0 - kp * kp );
			exp1 = exp(II*kz*dl);
			exp2 = exp(2.*II*kz*dl);
			reflTrans( eps[k], eps[k-1], kp/k0, pol, r, t );
			
			//printLayerInfo(k, eps[k], pol, r, t, exp1, exp2, kp, kz);

			s[ind11(k-1,l)] = s[ind11(k,l)] * exp1 * (1. - r*r) /( t*(1.-r*s[ind21(k,l)]) );
			s[ind12(k-1,l)] = s[ind12(k,l)] + r * s[ind11(k,l)] * s[ind22(k,l)] 
							/ ( 1. - r*s[ind21(k,l)] );
			s[ind21(k-1,l)] = exp2 * (s[ind21(k,l)] - r) / (1. - r*s[ind21(k,l)]);
			s[ind22(k-1,l)] = t * exp1 * s[ind22(k,l)] / ( 1. - r*s[ind21(k,l)] );
		}
	}
}

cdouble SMatrix::S11(int a, int b, int TM) {
	if(TM==1)
		return sTM[ind11(a,b)];
	else
		return sTE[ind11(a,b)];
}

cdouble SMatrix::S12(int a, int b, int TM) {
	if(TM==1)
		return sTM[ind12(a,b)];
	else
		return sTE[ind12(a,b)];
}

cdouble SMatrix::S21(int a, int b, int TM) {
	if(TM==1)
		return sTM[ind21(a,b)];
	else
		return sTE[ind21(a,b)];
}

cdouble SMatrix::S22(int a, int b, int TM) {
	if(TM==1)
		return sTM[ind22(a,b)];
	else
		return sTE[ind22(a,b)];
}

/*
int main() {
	// test Fresnel coefficients
	cdouble rs, rp, ts, tp;
	cdouble *epsA;
	cdouble eps0 = 1.0;
	cdouble eps1 = 3.5*3.5;
	reflTrans(eps0, eps1, sin(M_PI/6.), 0, rs, ts);
	reflTrans(eps0, eps1, sin(M_PI/6.), 1, rp, tp);
	std::cout << "rs: " << rs << "\nrp: " << rp << "\nts: " 
		<< ts << "\ntp: " << tp << "\n" << std::endl;

	// test scattering matrix
	double w = 2*M_PI*3e8/400e-9;
	int noLayers = 5;
	epsfn epsSi = [=] (double d) { return epsConst(d,cdouble(30.875,4.311)); };
	epsfn eps[] = {epsVac, epsSi, epsVac, epsSi, epsVac};
	epsA = new cdouble[noLayers];
	for(int i=0; i<noLayers; ++i)
		epsA[i] = eps[i](w);
	double d[] = { 0.25e-6, 0.5e-6, 0.25e-6 };
	SMatrix *S = new SMatrix(epsA, d, noLayers-1,w/3.e8,w/3.e8*sqrt(21./25.));

	int ind0, ind1;
	for(int pol=0; pol<=1; ++pol) {
		for(int ind=0; ind<=1; ++ind) {
			for(int i=0; i<noLayers; ++i) {
				ind0 = (ind==0) ? 0 : i;
				ind1 = (ind==0) ? i : noLayers-1;
				std::cout << "layer: " << i << "   pol: " << pol << std::endl;
				std::cout << "  ind0: " << ind0 << "  ind1: " << ind1 << std::endl;
				std::cout << "  intInd11: " << ind11(ind0,ind1) << "  intInd22: " << ind22(ind0,ind1) << std::endl;;
				std::cout << "  s11(ind0,ind1): " << S->S11(ind0,ind1,pol) << std::endl;
				std::cout << "  s12(ind0,ind1): " << S->S12(ind0,ind1,pol) << std::endl;
				std::cout << "  s21(ind0,ind1): " << S->S21(ind0,ind1,pol) << std::endl;
				std::cout << "  s22(ind0,ind1): " << S->S22(ind0,ind1,pol) << std::endl;
			}
		}
	}
}
*/
void eyeS(cdouble *S, int i) { // insert identity matrix
	S[i+0] = 1; S[i+1] = 0;
	S[i+2] = 0; S[i+3] = 1;
}
int ind11(int a, int b) {
	return 4*(a+b);
}
int ind12(int a, int b) {
	return 4*(a+b)+1;
}
int ind21(int a, int b) {
	return 4*(a+b)+2;
}
int ind22(int a, int b) {
	return 4*(a+b)+3;
}
