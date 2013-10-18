
#include <iostream>
#include <complex>
#include "SMatrix.hpp"
#include "multilayer.hpp"

typedef std::complex<double> cdouble;
const cdouble II(0.,1.);

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

// flux from zs to zl for given geometry (g) and kp, k0 (S)
//   nHat = 1. or -1.
double flux(const mlgeo *g, int l, int s, double zl, double zs, 
	double k0, double kp, double nHat) {
	if (g->eps[s].imag() == 0)
		return 0;

	pwaves pTE, pTM;
	cdouble hgf;

	SMatrix *S = new SMatrix(g, k0, kp);
	pWavesL(S, l, s, TE, pTE); 	
	pWavesL(S, l, s, TM, pTM);
	hgf = gfFlux(S, g, pTE, pTM, l, s, zl, zs);
	delete S;
	return nHat * std::real( S->k0 * S->k0 * II * g->eps[s].imag() * hgf / (M_PI * M_PI) );
}

// Green's functions contribution to flux
cdouble gfFlux(const SMatrix *S, const mlgeo *g, const pwaves &pTE, const pwaves &pTM, 
	int l, int s, double zl, double zs) {
			
	cdouble xl, xs, Ae, Am, Be, Bm, Ce, Cm, De, Dm;
	cdouble gEpp, gEpz, gEtt, gHpt, gHtp, gHtz; // gEzp, gEzz, gHzt

	const cdouble &kzl = S->kz[l];
	const cdouble &kzs = S->kz[s];
	const cdouble &kl = S->k[l];
	const cdouble &ks = S->k[s];
	const cdouble &kp = S->kp;
	xl = II * kzl * zl; 
	xs = II * kzs * (zs + g->z[s]); // abs, not rel pos.

	Ae = pTE.Al * exp(xl-xs);
	Be = pTE.Bl * exp(-xl-xs);
	Ce = pTE.Cl * exp(xl+xs);
	De = pTE.Dl * exp(xl-xs);
	Am = pTM.Al * exp(xl-xs);
	Bm = pTM.Bl * exp(-xl-xs);
	Cm = pTM.Cl * exp(xl+xs);
	Dm = pTM.Dl * exp(xl-xs);
	
	// electric dyadic Green's functions
	// NOTE: extra factor of kp (elsewhere in Francoeur Eqn.) to make g's unitless
	gEpp = II * kzl * kp / (2. * ks * kl) * ( Am - Bm - Cm + Dm );
	gEpz = II * kzl * kp * kp / (2. * kzs * ks * kl) * ( -Am + Bm - Cm + Dm ); 
	gEtt = II * kp / (2. * kzs) * ( Ae + Be + Ce + De );
	//gEzp = II * kp * kp / (2. * ks * kl) * ( -Am - Bm + Cm + Dm );
	//gEzz = II * kp * kp * kp / (2. * ksz * ks * kl) * ( Am + Bm + Cm + Dm );

	// magnetic dyadic Green's functions
	gHpt = kzl / (2. * kzs) * ( Ae - Be + Ce - De );
	gHtp = kl / (2. * ks) * ( -Am - Bm + Cm + Dm ); 
	gHtz = kl * kp / (2. * ks * kzs) * ( Am + Bm + Cm + Dm );
	//gHzt = kp / (2. * kzs) * ( -Ae - Be - Ce - De );

	return gEpp*std::conj(gHtp) + gEpz*std::conj(gHtz) - gEtt*std::conj(gHpt);
}

