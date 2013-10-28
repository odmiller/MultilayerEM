
#include <iostream>
#include <complex>
#include "SMatrix.hpp"
#include "mlgeo.hpp"
#include "fieldFunctions.hpp"

typedef std::complex<double> cdouble;
const cdouble II(0.,1.);
using std::conj;
using std::real;
using std::imag;

// Flux from s to zl for given geometry (g) and kp, k0
// Integrated over emitter evenly distributed in s
//   if don't want integral, gfFlux's below provide simple interface
//   nHat = 1. or -1.
//   in Francoeur et. al., zl --> zc and zs --> zprime
//   units of kp/k0 matter only relative to zl/zs
//   zl, zs are RELATIVE positions, starting from left boundary of resp. layers
//   not yet permitted: zl<0 if l==0, zl>0 if l==N (bndry condtions change)
//
// NOTE: Do not use this version of the function if you have more than one 
// emitting layer.  This version recomputes the scattering matrix every time.
// Precompute the scattering matrix (for a given kp & k0), then use the
// function below of the same name but different signature.
double flux(const mlgeo *g, double k0, double kp, int l, int s, 
			double zl, double nHat) {
	if ( imag(g->eps(s)) == 0)
		return 0;
	SMatrix *S = new SMatrix(g, k0, kp);
	double f = flux(S, l, s, zl, g->eps(s), g->d(s), nHat);
	delete S;
	return f;
}
double flux(const SMatrix *S, int l, int s, double zl,
			cdouble epsS, double ds, double nHat) {
	if ( imag(epsS) == 0)
		return 0;
	pwaves pTE, pTM;
	pWavesL(S, l, s, TE, pTE); 	
	pWavesL(S, l, s, TM, pTM);
	cdouble hgf = gfFluxTE(S, pTE, l, s, zl, ds, true); 
		 		+ gfFluxTM(S, pTM, l, s, zl, ds, true);  
	return nHat * real( S->k0 * S->k0 * II * imag(epsS) * hgf / (M_PI * M_PI) );
}

double meanEnergy(double w, double T) {
	const double hbar = 6.62606957e-34/2./M_PI;
	const double kB = 1.3806488e-23;
	return hbar*w / (exp(hbar*w/(kB*T)) - 1); // Boltzmann "Theta"
}

// compute just the partial waves in layer l
//   Splus and Sminus taken out, inserted in flux eqn. (as in Francoeur)
void pWavesL(const SMatrix *S, int l, int s, int pol, pwaves &p) {
	int N = S->N;
	// partial waves in s,0,N layers (A0=BN=C0=DN=0)
	cdouble Bs = S->S21(s,N,pol) / (1. - S->S12(0,s,pol)*S->S21(s,N,pol)) ; 
	cdouble As = S->S12(0,s,pol) * Bs;
	cdouble B0 = S->S22(0,s,pol) * Bs;
	cdouble AN = S->S11(s,N,pol) * (As + 1.);
	cdouble Cs = S->S12(0,s,pol) / (1. - S->S12(0,s,pol)*S->S21(s,N,pol));
	cdouble Ds = S->S21(s,N,pol) * Cs;
	cdouble CN = S->S11(s,N,pol) * Cs;
	cdouble D0 = S->S22(0,s,pol) * (Ds + 1.);
	if(l<s) { // flux to the left of emitter
		p.Bl = B0 / S->S22(0,l,pol);
		p.Al = S->S12(0,l,pol) * p.Bl;
		p.Dl = D0 / S->S22(0,l,pol);
		p.Cl = S->S12(0,l,pol) * p.Dl;
	} else { // flux to the right of emitter
		p.Al = AN / S->S11(l,N,pol);
		p.Bl = S->S21(l,N,pol) * p.Al;
		p.Cl = CN / S->S11(l,N,pol);
		p.Dl = S->S21(l,N,pol) * p.Cl;
	}
}

// integral of product of exponential terms arising in Green's function
// computations (i.e. integral of term in zsProd() from 0 to ds)
cdouble zsInt(int s1, int s2, cdouble kzs, double ds) {
	if(s1==1 && s2==1)
		return (exp(2. * II * real(kzs) * ds) - 1.) / (2. * II * real(kzs));
	else if(s1==1 && s2==-1)
		return (1. - exp(-2. * imag(kzs) * ds)) / (2. * imag(kzs));
	else if(s1==-1 && s2==1)
		return (exp(2. * imag(kzs) * ds) - 1.) / (2. * imag(kzs));
	else if(s1==-1 && s2==-1)
		return (1. - exp(-2. * II * real(kzs) * ds)) / (2. * II * real(kzs));
	return -1; // shouldn't get here
}

// product of exponential terms arising in Green's function computations
cdouble zsProd(int s1, int s2, cdouble kzs, double zs) {
	return exp(II * zs * (double(s1) * kzs + double(s2) * conj(kzs)));
}

// TE flux from Green's functions
// if integrate==true (default), then integral over emitter layer 
//   done analytically. In this case xs = thickness of layer s
// if integrate==false, then xs is the location of the emitter
cdouble gfFluxTE(const SMatrix *S, const pwaves &pTE, 
				int l, int s, double zl, double xs, bool integrate) {
	cdouble kzl = S->kz[l];
	cdouble kzs = S->kz[s];
	cdouble kp = S->kp;
	cdouble xl = II * kzl * zl; 
	cdouble A = pTE.Al * exp(xl);
	cdouble B = pTE.Bl * exp(-xl);
	cdouble C = pTE.Cl * exp(xl);
	cdouble D = pTE.Dl * exp(-xl);

	cdouble fTE = 0;
	int gES[4] = {-1,-1,+1,+1}; // signs in exponent terms
	int gHS[4] = {+1,+1,-1,-1};

	cdouble (*spaceFx)(int, int, cdouble, double);
	if (integrate)
		spaceFx = zsInt;
	else
		spaceFx = zsProd;

	// Note the overall neg. sign below
	cdouble prefac = II * kp / (4. * kzs) * conj(kzl / kzs);
	cdouble gEtt[4] = {A, B, C, D};
	cdouble gHpt[4] = {A, -B, C, -D};
	for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			fTE -= prefac * gEtt[i] * conj(gHpt[j]) * spaceFx(gES[i], gHS[j], kzs, xs);
	
	return fTE;
}

// TM flux from Green's functions
cdouble gfFluxTM(const SMatrix *S, const pwaves &pTM, 
				int l, int s, double zl, double xs, bool integrate) {
	cdouble kzl = S->kz[l];
	cdouble kzs = S->kz[s];
	cdouble kl = S->k[l];
	cdouble ks = S->k[s];
	cdouble kp = S->kp;
	cdouble xl = II * kzl * zl; 
	cdouble A = pTM.Al * exp(xl);
	cdouble B = pTM.Bl * exp(-xl);
	cdouble C = pTM.Cl * exp(xl);
	cdouble D = pTM.Dl * exp(-xl);

	cdouble fTM = 0;
	int gES[4] = {-1,-1,+1,+1}; // signs in exponent terms
	int gHS[4] = {+1,+1,-1,-1};

	cdouble (*spaceFx)(int, int, cdouble, double);
	if (integrate)
		spaceFx = zsInt;
	else
		spaceFx = zsProd;

	// TM1
	cdouble prefac = II * kzl * kp / (4. * ks * kl) * conj(kl / ks);
	cdouble gEpp[4] = {A, -B, -C, D};
	cdouble gHtp[4] = {-A, -B, C, D};
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			fTM += prefac * gEpp[i] * conj(gHtp[j]) * spaceFx(gES[i], gHS[j], kzs, xs); 

	// TM2
	prefac *= kp * conj(kp) / (kzs * conj(kzs));
	cdouble gEpz[4] = {-A, B, -C, D};
	cdouble gHtz[4] = {A, B, C, D};
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			fTM += prefac * gEpz[i] * conj(gHtz[j]) * spaceFx(gES[i], gHS[j], kzs, xs);
	
	return fTM;
}

// DELETE THIS ONCE IT HAS BEEN TESTED AGAINST!!
// Green's functions contribution to flux, from single zs in s
// Very important: zprime in Francoeur et. al. should be relative 
// to zs (Francoeur notation), not the absolute position.
cdouble gfFluxSP(const SMatrix *S, const pwaves &pTE, const pwaves &pTM, 
	int l, int s, double zl, double zs) {
			
	cdouble xl, xs, Ae, Am, Be, Bm, Ce, Cm, De, Dm;
	cdouble gEpp, gEpz, gEtt, gHpt, gHtp, gHtz; // gEzp, gEzz, gHzt

	const cdouble &kzl = S->kz[l];
	const cdouble &kzs = S->kz[s];
	const cdouble &kl = S->k[l];
	const cdouble &ks = S->k[s];
	const cdouble &kp = S->kp;
	xl = II * kzl * zl; 
	xs = II * kzs * zs; 

	Ae = pTE.Al * exp(xl-xs);
	Be = pTE.Bl * exp(-xl-xs);
	Ce = pTE.Cl * exp(xl+xs);
	De = pTE.Dl * exp(-xl+xs);
	Am = pTM.Al * exp(xl-xs);
	Bm = pTM.Bl * exp(-xl-xs);
	Cm = pTM.Cl * exp(xl+xs);
	Dm = pTM.Dl * exp(-xl+xs);
	
	// electric dyadic Green's functions
	// NOTE: extra factor of kp (elsewhere in Francoeur Eqn.) to make g's unitless
	gEpp = II * kzl * kp / (2. * ks * kl) * ( Am - Bm - Cm + Dm );
	gEpz = II * kzl * kp * kp / (2. * kzs * ks * kl) * ( -Am + Bm - Cm + Dm ); 
	gEtt = II * kp / (2. * kzs) * ( Ae + Be + Ce + De );

	// magnetic dyadic Green's functions
	gHpt = kzl / (2. * kzs) * ( Ae - Be + Ce - De );
	gHtp = kl / (2. * ks) * ( -Am - Bm + Cm + Dm ); 
	gHtz = kl * kp / (2. * ks * kzs) * ( Am + Bm + Cm + Dm );

	return gEpp*conj(gHtp) + gEpz*conj(gHtz) - gEtt*conj(gHpt);
}

