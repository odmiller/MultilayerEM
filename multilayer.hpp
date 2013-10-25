
#include <complex>
#include "SMatrix.hpp"

// primary function: flux from emitter at zs to zl (in layers s/l resp.)
//   at free-space wavevector k0 and parallel wavevector k0
// zs<0 indicates integration over emitter layer
// output has units: number / distance^2 / wavevector / freq
double flux(const mlgeo *g, int l, int s, double zl,
	double k0, double kp, double nHat, double zs=-1.);

// mean energy (units = J) of Planck oscillator at freq w, temp. T
double meanEnergy(double w, double T);

// internal struct's / functions used 
//   to compute Green's functions & partial waves
struct pwaves {
	std::complex<double> Al, Bl, Cl, Dl;
};

// emission from all of s, analytically integrated
std::complex<double> gfFlux(const SMatrix *S, const pwaves &pTE,
	const pwaves &pTM, int l, int s, double zl, double ts);

// function to compute emission from a single point (i.e. not integrated
// over whole emitter layer).  Usually this function should not be needed.
// Note that the signatures are identical.
std::complex<double> gfFluxSP(const SMatrix *S, const pwaves &pTE,
	const pwaves &pTM, int l, int s, double zl, double zs);

// compute partial waves l relative to a source in s
void pWavesL(const SMatrix *S, int l, int s, int pol, pwaves &p);
