
#include <iostream>
#include <complex>
#include "SMatrix.hpp"
#include "mlgeo.hpp"
#include "multilayer.hpp"

typedef std::complex<double> cdouble;
const cdouble II(0.,1.);
using std::conj;
using std::real;
using std::imag;

// for given stack (sMatrix:S,kz,k), compute flux at zl (in layer l) from 
// emitter at zs (in layer s, with permittivity epsS), for k0 and kp (in S)
//   in Francoeur et. al., zl --> zc and zs --> zprime
// units of kp,k0 matter only relative to zl,zs
// zl, zs are RELATIVE positions, starting from left boundary of resp. layers
//   zl \in [0,t[l]] if 0<l<N
//  not yet permitted: zl<0 if l==0, zl>0 if l==N (bndry condtions change)

// compute just the partial waves in layer l
//   should have update() function for S-Matrix
//   Splus and Sminus taken out, inserted in flux eqn. (as in Ref.)
void pWavesL(const SMatrix *S, int l, int s, int pol, pwaves &p) {
	cdouble As, AN, B0, Bs, Cs, CN, D0, Ds;
	int N = S->N;

	// partial waves in s,0,N layers (A0=BN=C0=DN=0)
	Bs = S->S21(s,N,pol) / (1. - S->S12(0,s,pol)*S->S21(s,N,pol)) ; 
	As = S->S12(0,s,pol) * Bs;
	B0 = S->S22(0,s,pol) * Bs;
	AN = S->S11(s,N,pol) * (As + 1.);
	Cs = S->S12(0,s,pol) / (1. - S->S12(0,s,pol)*S->S21(s,N,pol));
	Ds = S->S21(s,N,pol) * Cs;
	CN = S->S11(s,N,pol) * Cs;
	D0 = S->S22(0,s,pol) * (Ds + 1.);
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

// flux from s to zl for given geometry (g) and kp, k0 (S)
// either from a single zs, or from entire s layer
//   nHat = 1. or -1.
double flux(const mlgeo *g, int l, int s, double zl,  
	double k0, double kp, double nHat, double zs) {
	if ( imag(g->eps(s)) == 0)
		return 0;

	pwaves pTE, pTM;
	cdouble hgf;

	SMatrix *S = new SMatrix(g, k0, kp);
	pWavesL(S, l, s, TE, pTE); 	
	pWavesL(S, l, s, TM, pTM);

	if(zs<0) // not allowed otherwise (use to signify layer calc - default val=-1)
		hgf = gfFlux(S, pTE, pTM, l, s, zl, g->d(s));
	else
		hgf = gfFluxSP(S, pTE, pTM, l, s, zl, zs);
	delete S;
	return nHat * real( k0 * k0 * II * imag(g->eps(s)) * hgf / (M_PI * M_PI) );
}

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

// Should split this up into TE/TM (also merge with non-integrated version)
// Green's functions contributions to flux, entire s layers
// integral over emitter layer done analytically
cdouble gfFlux(const SMatrix *S, const pwaves &pTE, const pwaves &pTM, 
	int l, int s, double zl, double ts) {

	cdouble xl, xs, Ae, Am, Be, Bm, Ce, Cm, De, Dm;
	cdouble gEpp[4], gEpz[4], gEtt[4], gHpt[4], gHtp[4], gHtz[4]; // gEzp, gEzz, gHzt

	const cdouble &kzl = S->kz[l];
	const cdouble &kzs = S->kz[s];
	const cdouble &kl = S->k[l];
	const cdouble &ks = S->k[s];
	const cdouble &kp = S->kp;
	xl = II * kzl * zl; 

	Ae = pTE.Al * exp(xl);
	Be = pTE.Bl * exp(-xl);
	Ce = pTE.Cl * exp(xl);
	De = pTE.Dl * exp(-xl);
	Am = pTM.Al * exp(xl);
	Bm = pTM.Bl * exp(-xl);
	Cm = pTM.Cl * exp(xl);
	Dm = pTM.Dl * exp(-xl);

	cdouble fTM1 = 0, fTM2 = 0, fTE = 0;
	int gES[4] = {-1,-1,+1,+1};
	int gHS[4] = {+1,+1,-1,-1};

	// TM (1)
	//std::cout << "kzl*kp/ks*ks " << kzl*kp/(ks*ks) << std::endl;
	gEpp[0] = II * kzl * kp / (2. * ks * kl) * Am;
	gEpp[1] = II * kzl * kp / (2. * ks * kl) * (-Bm);
	gEpp[2] = II * kzl * kp / (2. * ks * kl) * (-Cm);
	gEpp[3] = II * kzl * kp / (2. * ks * kl) * Dm;
	gHtp[0] = kl / (2. * ks) * (-Am);
	gHtp[1] = kl / (2. * ks) * (-Bm);
	gHtp[2] = kl / (2. * ks) * Cm;
	gHtp[3] = kl / (2. * ks) * Dm;
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			fTM1 += gEpp[i] * conj(gHtp[j]) * zsInt(gES[i], gHS[j], kzs, ts); 

	// TM (2)
	gEpz[0] = II * kzl * kp * kp / (2. * kzs * ks * kl) * (-Am);
	gEpz[1] = II * kzl * kp * kp / (2. * kzs * ks * kl) * Bm;
	gEpz[2] = II * kzl * kp * kp / (2. * kzs * ks * kl) * (-Cm);
	gEpz[3] = II * kzl * kp * kp / (2. * kzs * ks * kl) * Dm;
	gHtz[0] = kl * kp / (2. * ks * kzs) * Am;
	gHtz[1] = kl * kp / (2. * ks * kzs) * Bm;
	gHtz[2] = kl * kp / (2. * ks * kzs) * Cm;
	gHtz[3] = kl * kp / (2. * ks * kzs) * Dm;
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			fTM2 += gEpz[i] * conj(gHtz[j]) * zsInt(gES[i], gHS[j], kzs, ts);

	// TE pol: note the overall neg. sign below
	gEtt[0] = II * kp / (2. * kzs) * Ae;
	gEtt[1] = II * kp / (2. * kzs) * Be;
	gEtt[2] = II * kp / (2. * kzs) * Ce;
	gEtt[3] = II * kp / (2. * kzs) * De;
	gHpt[0] = kzl / (2. * kzs) * Ae;
	gHpt[1] = kzl / (2. * kzs) * (-Be);
	gHpt[2] = kzl / (2. * kzs) * Ce;
	gHpt[3] = kzl / (2. * kzs) * (-De);
	for(int i=0; i<4; ++i)
		for(int j=0; j<4; ++j)
			fTE -= gEtt[i] * conj(gHpt[j]) * zsInt(gES[i], gHS[j], kzs, ts);

	return fTM1+fTM2+fTE;
}

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

