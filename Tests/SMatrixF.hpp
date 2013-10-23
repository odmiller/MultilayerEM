
#ifndef SMATRIXF_H
#define SMATRIXF_H

#include <complex>
#include "SMatrix.hpp"

// SMatrixF contains the scattering matrices (and indexing functions)
//   for a multilayer stack with permittivity eps and thicknesses t
class SMatrixF {
	public:
		SMatrixF(const std::complex<double> *eps, const double *d, 
				int NIn, double k0In, double kpIn);
		SMatrixF(const mlgeo *g, double k0In, double kpIn);
		~SMatrixF();
		std::complex<double> S11(int a, int b, int pol) const; 
		std::complex<double> S12(int a, int b, int pol) const;
		std::complex<double> S21(int a, int b, int pol) const;
		std::complex<double> S22(int a, int b, int pol) const;
		void printSMatrix();
		std::complex<double> *kz; // kz = sqrt( eps*k0^2 - kp^2 )
		std::complex<double> *k; // sqrt(eps)*k0	
		double kp, k0; // par wavevector, k0
		int N;

	private:
		void initMatrix(const std::complex<double> *eps, const double *d);
		std::complex<double> *sTE;
		std::complex<double> *sTM;
};

#endif
