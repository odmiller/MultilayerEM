
// MULTI-LAYER GEOMETRY CLASS

#ifndef MLGEO_H
#define MLGEO_H

#include <complex>

class mlgeo {
	public:
		mlgeo(const std::complex<double> *epsIn,const double *dIn, int NIn);
		~mlgeo();
		std::complex<double> eps(int i) const;
		double d(int i) const;
		double z(int i) const;
		void print() const;
		int N;
	private:
		std::complex<double> *eps_;
		double *d_;
		double *z_;
};

#endif
