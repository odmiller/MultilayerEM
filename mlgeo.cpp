
#include <assert.h>
#include <iostream>
#include <complex>
#include "mlgeo.hpp"

typedef std::complex<double> cdouble;

mlgeo::mlgeo(const cdouble *epsIn, const double *dIn, int NIn) 
	: N(NIn), eps_(new cdouble[NIn+1]), d_(new double[NIn-1]), z_(new double[NIn]) {
	z_[0] = 0.; // set first interface at zero
	for(int i=0; i<=N; ++i) {
		eps_[i] = epsIn[i];
		if(i<N-1)
			d_[i] = dIn[i];
		if(i>0 && i<N)
			z_[i] = z_[i-1] + d_[i-1];
	}
}

mlgeo::~mlgeo() {
	delete[] eps_;
	delete[] d_;
	delete[] z_;
}

cdouble mlgeo::eps(int i) const {
	return eps_[i];
}

double mlgeo::d(int i) const {
	assert(i>0 && i<N);
	return d_[i-1]; 
}

double mlgeo::z(int i) const {
	assert(i>0);
	return z_[i-1];
}

void mlgeo::print() const {
	std::cout << "-------------------------" << std::endl;
	std::cout << "Number of interfaces: " << N << std::endl;
	std::cout << "Number of layers: " << N + 1 << std::endl;
	for(int i=0; i<=N; ++i) {
		std::cout << "layer: " << i << std::endl;
		std::cout << "  eps: " << eps(i) << std::endl;
		if(i>0 && i<N)
			std::cout << "  d: " << d(i) << std::endl;
		if(i>0)
			std::cout << "  z: " << z(i) << std::endl;
	}
	std::cout << "-------------------------" << std::endl;
}


