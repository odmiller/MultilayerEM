
#ifndef SMATRIX_H
#define SMATRIX_H

#include <complex>

// SMatrix contains the scattering matrices (and indexing functions)
//   for a multilayer stack with permittivity eps and thicknesses t
class SMatrix {
	public:
		SMatrix(std::complex<double> *eps, double *d, int N, double k0, double kp);
		std::complex<double> S11(int a, int b, int TM);
		std::complex<double> S12(int a, int b, int TM);
		std::complex<double> S21(int a, int b, int TM);
		std::complex<double> S22(int a, int b, int TM);
/*		std::complex<double> *kz; // kz = sqrt( eps*k0^2 - kp^2 )
		std::complex<double> *k; // sqrt(eps)*k0	
		std::complex<double> kp; // par wavevector*/

	private:
		std::complex<double> *sTE;
		std::complex<double> *sTM;
};

#endif
