
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

const double hbar = 6.62606957e-34/2./M_PI;
const double kB = 1.3806488e-23;
const double c0 = 299792458.;
const double Z0 = 119.9169832 * M_PI;

// Flux from s to zl for given geometry (g) and kp, k0
// Integrated over emitter evenly distributed in s
//   if don't want integral, gfFlux's below provide simple interface
//   nHat = 1. or -1.
//   in Francoeur et. al., zl --> zc and zs --> zprime
//   units of kp/k0 matter only relative to zl/zs
//   zl, zs are RELATIVE positions, starting from left boundary of resp. layers
//   not yet permitted: zl<0 if l==0, zl>0 if l==N (bndry condtions change)
//
// Specify all emitting layers in s array.  First version computes/discards the 
//   scattering matrix for you, second accepts precomputed S.  nHat is assumed 
//   to be the same for each emittering layer
double flux(const mlgeo &g, double k0, double kp, int l, double zl, 
			const int *s, int Ns, double nHat) { 
	SMatrix S = SMatrix(g, k0, kp);
	return flux(g, S, l, zl, s, Ns, nHat);
}
double flux(const mlgeo &g, const SMatrix &S, int l, double zl, 
			const int *s, int Ns, double nHat) { 
	pwaves pTE, pTM;
	double f = 0;
	for (int sind = 0, si = s[sind]; sind < Ns; ++sind, ++si) {
		if (imag(g.eps(si)) == 0)
			continue;
		pWavesL(S, l, si, TE, &pTE);
		pWavesL(S, l, si, TM, &pTM);
		f -= nHat * S.k0 * S.k0 / (M_PI * M_PI) * imag(g.eps(si)) 
			* imag( gfFluxTE(S, pTE, l, si, zl, g.d(si), nHat)
			+ gfFluxTM(S, pTM, l, si, zl, g.d(si), nHat) );
	}
	return f;
}

// flux of a blackbody (/ dist^2 / freq), not / wavevector!
// multiply by (1/a)^2 to get 1/m^2
double flux_bb(double k0) {
	return k0 * k0 / (4. * M_PI * M_PI);
}

// Integrated flux for blackbody (in vacuum)
// Given by Stefan-Boltzmann (or integral of mean_energy * dos_bb)
// multiply by hbar*c^2/lscale^4 to get W/m^2
double flux_bb_int(double lscale, double T) {
	return pow(lscale * kB * T / (hbar * c0), 4) * M_PI * M_PI / 60.;
}

// density of states, DOS (technically electric DOS, i.e. DOS of an electric dipole)
//   cf. e.g. Joulain et. al. PRB 68, 245405 (2003)
double dos(const mlgeo &g, double k0, double kp, int l, double zl) {
	SMatrix S = SMatrix(g, k0, kp);
	return dos(S, l, zl);
}

// same as above but with S precomputed (e.g. if multiple zl's desired)
double dos(const SMatrix &S, int l, double zl) {
	pwaves pTE, pTM;	
	pWavesL(S, l, l, TE, &pTE); // source layer = emitter layer
	pWavesL(S, l, l, TM, &pTM);

	cdouble kzl = S.kz[l]; // kzs = kzl
	cdouble kl = S.k[l];
	cdouble kp = S.kp;
	cdouble xl = II * kzl * zl; // xs = xl
	
	cdouble Ae = pTE.Al;
	cdouble Be = pTE.Bl * exp(-2.*xl);
	cdouble Ce = pTE.Cl * exp(2.*xl);
	cdouble De = pTE.Dl;
	cdouble Am = pTM.Al;
	cdouble Bm = pTM.Bl * exp(-2.*xl);
	cdouble Cm = pTM.Cl * exp(2.*xl);
	cdouble Dm = pTM.Dl;
		
	// extra factors of 1 are extra source term when l==s
	cdouble Epp = II * kzl * kp / (kl * kl) 
				* (Am - Bm - Cm + Dm + 1.); 
	cdouble Ett = II * kp / kzl
				* (Ae + Be + Ce + De + 1.);
	cdouble Ezz = II * kp * kp * kp / (kzl * kl * kl) 
				* (Am + Bm + Cm + Dm + 1.);
	return S.k0 * S.k0 * imag(Epp + Ett + Ezz) / (2. * M_PI * M_PI);   
}

double dos_vacuum(double k0) {
	return k0 * k0 / (2. * M_PI * M_PI);
}

// multiply by hbar*c/lscale to get J
double mean_energy(double k0, double lscale, double T) {
	return k0 / (exp(hbar * c0 * k0 / (lscale * kB * T)) - 1); // Boltzmann "Theta"
}

void reflTrans(const mlgeo &g, double k0, double theta, int pol, 
				cdouble *r, double *R, cdouble *t, double *T) {
	SMatrix S = new SMatrix(g, k0, k0 * sin(theta));
	*r = S->S21(0, g.N, pol);
	*t = S->S11(0, g.N, pol);
	*R = (*r) * conj(*r); 
	*T = real(sqrt(g.eps(N))) / real(sqrt(g.eps(0))) * (*t) * conj(*t); 
}

// compute just the partial waves in layer l
//   Splus and Sminus taken out, inserted in flux eqn. (as in Francoeur)
void pWavesL(const SMatrix &S, int l, int s, int pol, pwaves *p) {
	int N = S.N;
	if (s==0) { // emitting half-space
		p->Cl = 0; // no waves emitted in the backward direction
		p->Dl = 0;
		if (l==0) {
			p->Al = 0;
			p->Bl = S.S21(l,N,pol);
		} else {
			p->Al = S.S11(0,N,pol) / S.S11(l,N,pol);
			p->Bl = S.S21(l,N,pol) * p->Al;
		}
	} else if (s==N) {
		p->Al = 0; // no waves in forward dir
		p->Bl = 0;
		if (l==N) {
			p->Dl = 0;
			p->Cl = S.S12(0,l,pol);
		} else {
			p->Dl = S.S22(0,N,pol) / S.S22(0,l,pol);
			p->Cl = S.S12(0,l,pol) * p->Dl;
		}
	} else {
		// partial waves in s,0,N layers (A0=BN=C0=DN=0)
		cdouble Bs = S.S21(s,N,pol) / (1. - S.S12(0,s,pol)*S.S21(s,N,pol)) ; 
		cdouble As = S.S12(0,s,pol) * Bs;
		cdouble B0 = S.S22(0,s,pol) * Bs;
		cdouble AN = S.S11(s,N,pol) * (As + 1.);
		cdouble Cs = S.S12(0,s,pol) / (1. - S.S12(0,s,pol)*S.S21(s,N,pol));
		cdouble Ds = S.S21(s,N,pol) * Cs;
		cdouble CN = S.S11(s,N,pol) * Cs;
		cdouble D0 = S.S22(0,s,pol) * (Ds + 1.);
		if (l==s) { // flux, emitter in same layer
			p->Al = As;
			p->Bl = Bs;
			p->Cl = Cs;
			p->Dl = Ds;
		} else if (l<s) { // flux to the left of emitter
			p->Bl = B0 / S.S22(0,l,pol);
			p->Al = S.S12(0,l,pol) * p->Bl;
			p->Dl = D0 / S.S22(0,l,pol);
			p->Cl = S.S12(0,l,pol) * p->Dl;
		} else { // flux to the right of emitter
			p->Al = AN / S.S11(l,N,pol);
			p->Bl = S.S21(l,N,pol) * p->Al;
			p->Cl = CN / S.S11(l,N,pol);
			p->Dl = S.S21(l,N,pol) * p->Cl;
		}
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
//   extra factor of kp (elsewhere in the integrand) to make dimensionless
// if integrate==true (default), then integral over emitter layer 
//   done analytically. In this case xs = thickness of layer s
// if integrate==false, then xs is the location of the emitter
cdouble gfFluxTE(const SMatrix &S, const pwaves &pTE, 
				int l, int s, double zl, double xs, bool integrate) {
	cdouble kzl = S.kz[l];
	cdouble kzs = S.kz[s];
	cdouble kp = S.kp;
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
cdouble gfFluxTM(const SMatrix &S, const pwaves &pTM, 
				int l, int s, double zl, double xs, bool integrate) {
	cdouble kzl = S.kz[l];
	cdouble kzs = S.kz[s];
	cdouble kl = S.k[l];
	cdouble ks = S.k[s];
	cdouble kp = S.kp;
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
