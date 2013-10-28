
#include <complex>
#include "SMatrix.hpp"

// primary function: flux from emitters in s to zl in layer l
//   at free-space wavevector k0 and parallel wavevector k0
// output has units: number / distance^2 / wavevector / freq
double flux(const mlgeo *g, double k0, double kp, int l, int s, 
	double zl, double nHat);

// same as above but with scattering matrix S precomputed
double flux(const SMatrix *S, int l, int s, double zl, 
	std::complex<double> epsS, double ds, double nHat);

// mean energy (units = Joules) of Planck oscillator at freq w, temp. T
double meanEnergy(double w, double T);


// Potentially useful struct's / functions:

struct pwaves {
	std::complex<double> Al, Bl, Cl, Dl;
};

// compute partial waves l relative to a source in layer s
void pWavesL(const SMatrix *S, int l, int s, int pol, pwaves &p);

// TE/TM flux rates at l from emitter in s (either all of s or a single point)
// Use these functions if don't want integral over emitter layer,
// or want polarization specificity
std::complex<double> gfFluxTE(const SMatrix *S, const pwaves &pTE,
	int l, int s, double zl, double xs, bool integrate=true);

std::complex<double> gfFluxTM(const SMatrix *S, const pwaves &pTM,
	int l, int s, double zl, double xs, bool integrate=true);

// function to compute emission from a single point (i.e. not integrated
// over whole emitter layer).  Usually this function should not be needed.
// Note that the signatures are identical.
std::complex<double> gfFluxSP(const SMatrix *S, const pwaves &pTE,
	const pwaves &pTM, int l, int s, double zl, double zs);

