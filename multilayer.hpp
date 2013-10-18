
#include <complex>
#include "SMatrix.hpp"

// primary function: flux from emitter at zs to zl (in layers s/l resp.)
//   at free-space wavevector k0 and parallel wavevector k0
double flux(const mlgeo *g, int l, int s, double zl, double zs, 
	double k0, double kp, double nHat);

// internal struct's / functions used 
//   to compute Green's functions & partial waves
struct pwaves {
	std::complex<double> Al, Bl, Cl, Dl;
};

std::complex<double> gfFlux(const SMatrix *S, const mlgeo *g, const pwaves &pTE,
	const pwaves &pTM, int l, int s, double zl, double zs);

void pWavesL(const SMatrix *S, int l, int s, int pol, pwaves &p);
