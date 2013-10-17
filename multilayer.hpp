
#include <complex>
#include "SMatrix.hpp"

struct pwaves {
	std::complex<double> Al, Bl, Cl, Dl;
}

void pWavesL(SMatrix *S, int l, int s, double zs, std::complex<double> epsS, 
	double kp, double k0, int pol, pwaves &p);
