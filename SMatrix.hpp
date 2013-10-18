
#ifndef SMATRIX_H
#define SMATRIX_H

#include <complex>

struct mlgeo {
	std::complex<double> *eps;
	double *d;
	double *z;
	int N;
};

const int TE = 1;
const int TM = 0;

// SMatrix contains the scattering matrices (and indexing functions)
//   for a multilayer stack with permittivity eps and thicknesses t
class SMatrix {
	public:
		SMatrix(const std::complex<double> *eps, const double *d, 
				int NIn, double k0In, double kpIn);
		SMatrix(const mlgeo *g, double k0In, double kpIn);
		~SMatrix();
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

void initGeo(const std::complex<double> *epsIn, const double *dIn, int NIn, mlgeo *g);
void printGeo(mlgeo *g);

// Fresnel coefficiencts at an interface
void reflTrans(std::complex<double> eps0, std::complex<double> eps1, 
	double kp, int pol, std::complex<double> &r, std::complex<double> &t);

#endif
