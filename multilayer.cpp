
#include <iostream>
#include <complex>
#include "SMatrix.hpp"

typedef std::complex<double> cdouble;
typedef cdouble II(0.,1.);

// for given stack (sMatrix), compute flux at zl (in layer l) from 
// emitter at zs (in layer s, with permittivity epsS), for angular 
// freq w and parallel wavevector kp
//   in Francoeur et. al., zl --> zc and zs --> zprime
// units of kp,k0 matter only relative to zl,zs
// zl, zs are positions relative to left boundary of resp. layers
//   zl \in [0,t[l]] if 0<l<N
//  not yet permitted: zl<0 if l==0, zl>0 if l==N (bndry condtions change)

// compute just the partial waves in layer l
void pWavesL(SMatrix *S, int l, int s, double zs, cdouble epsS, 
	double kp, double k0, int pol, pwaves &p) {
	
	cdouble kzs, Splus, Sminus;
	kzs = sqrt( epsS * k0 * k0 - kp * kp ); 
	Splus = exp(-1.*II*kzs*zs);
	Sminus = exp(II*kzs*zs);

	// partial waves in s,0,N layers (A0=BN=C0=DN=0)
	Bs = S->S21(s,N,pol) * Splus /( 1 - S->S12(0,s,pol)*S->S21(s,N,pol) ); 
	As = S->S12(0,s,pol) * Bs;
	B0 = S->S22(0,s,pol) * Bs;
	AN = S->S11(s,N,pol) * (As + Splus);
	Cs = S->S12(0,s,pol) * Sminus /( 1 - S->S12(0,s,pol)*S->S21(s,N,pol) );
	Ds = S->S21(s,N,pol) * Cs;
	CN = S->S11(s,N,pol) * Cs;
	D0 = S->S22(0,s,pol) * (Ds + Sminus);
	if(l<s) { // flux to the left of emitter
		p.Bl = B0 / S->S22(0,l,pol);
		p.Al = S->S12(0,l,pol) * Bl;
		p.Dl = D0 / S->S22(0,l,pol);
		p.Cl = S->S12(0,l,pol) * Dl;
	} else { // flux to the right of emitter
		p.Al = AN / S->S11(l,N,pol);
		p.Bl = S->S21(l,N,pol) * Al;
		p.Cl = CN / S->S11(l,N,pol);
		p.Dl = S->S21(l,N,pol) * Cl;
	}
}

// compute flux from partial waves
// TODO: check that imag(permittivity) > 0
double flux(pwaves &pTE, pwaves &pTM, double zl, double zs, double k0) {
			
	cdouble xl, xs, ks, kl;
	cdouble gEslpp, gEslpz, gEsltt;
	xl = sqrt( epsL * k0 * k0 - kp * kp ) * zl; // kzl*zl 
	xs = sqrt( epsS * k0 * k0 - kp * kp ) * zs; 


	// gEslpp (Francoeur notation)
	gEslpp = II * kzl / (2. * ks * kl) * (pTM.Al * 
}

cdouble gEslpp
