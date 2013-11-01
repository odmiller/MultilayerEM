
#include <complex>
#include "SMatrix.hpp"

/*****  primary functions  ********/

// flux from emitters in s to zl in layer l
//   at free-space wavevector k0 and parallel wavevector k0
// (number / distance^2 / wavevector / freq)
double flux(const mlgeo &g, double k0, double kp, int l, double zl,
	const int *s, int Ns, double nHat);
double flux(const mlgeo &g, const SMatrix &S, int l, double zl, 
	const int *s, int Ns, double nHat);

double flux_bb(double k0);
double flux_bb_int(double lscale, double T); // integrated over w at temp. T

// density of states (/ wavevector / distance^3 / second / freq)
double dos(const mlgeo &g, double k0, double kp, int l, double zl);
double dos(const SMatrix &S, int l, double zl);

double dos_vacuum(double k0); 

// mean energy of Planck oscillator at freq w, temp. T
double mean_energy(double k0, double lscale, double T);

// r, t = reflection / transmission coefficients
// R, T = reflected / transmitted power (divided by inc. power)
// at angle theta, for given polarization
void reflTrans(const mlgeo &g, double k0, double theta, int pol,
			   std::complex<double> *r, double *R, std::complex<double> *t, double *T);

// normalization for flux: 1/a^2 (i.e. a^2 * flux is in SI units)
// normalization for dos: 1/a^2/c
// 		where zl = zl_SI / a;
// mean_energy is already in J
//
/***** potentially useful structs / functions, innards ********/

struct pwaves {
	std::complex<double> Al, Bl, Cl, Dl;
};

// compute partial waves l relative to a source in layer s
void pWavesL(const SMatrix &S, int l, int s, int pol, pwaves *p);

// TE/TM flux rates at l from emitter in s (either all of s or a single point)
// Use these functions if don't want integral over emitter layer,
// or want polarization specificity
std::complex<double> gfFluxTE(const SMatrix &S, const pwaves &pTE,
	int l, int s, double zl, double xs, bool integrate=true);

std::complex<double> gfFluxTM(const SMatrix &S, const pwaves &pTM,
	int l, int s, double zl, double xs, bool integrate=true);
