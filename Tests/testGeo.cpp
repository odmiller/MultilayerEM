
#include <iostream>
#include <complex>
#include "mlgeo.hpp"
#include "SMatrix.hpp"

typedef std::complex<double> cdouble;

int main() {
	cdouble *eps = new cdouble[5];
	double *d = new double[3];
	
	eps[0] = 1;
	eps[1] = cdouble(-10.,1.);
	eps[2] = 1.;
	eps[3] = cdouble(-5.,2.);
	eps[4] = 1.;
	d[0] = 100;
	d[1] = 50;
	d[2] = 150;
	
	mlgeo *g = new mlgeo(eps, d, 4);
	g->print();
	
	double k0 = 2.*M_PI/400.;
	double kp = 0.;
	SMatrix *S = new SMatrix(g, k0, kp);	
	S->print();

	delete g;
	delete S;
	return 0;
}
