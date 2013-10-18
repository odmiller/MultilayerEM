
#include <iostream>
#include <complex>
#include "cubature.h"
#include "multilayer.hpp"
#include "SMatrix.hpp"

struct geoData {
	mlgeo *g;
	int l, s;
	double zl, k0, kp, nHat;
};

int intFx(unsigned ndim, const double *zs, void *fdata, unsigned fdim, double *fval) {
	geoData *f = (geoData *) fdata;
	*fval = flux(f->g, f->l, f->s, f->zl, *zs, f->k0, f->kp, f->nHat);
	return 0;
}

double fluxZInt(mlgeo *g, int l, int s, double zl, double k0, double kp, double nHat) {
	int res;
	unsigned fdim = 1, dim = 1;
	double xmin = 0, xmax = g->d[s-1], val, err; // s cannot equal 0 or N 
	double absError = 0, relError = 1e-4;
	size_t maxEval = 1e5;
	geoData fdata = {g, l, s, zl, k0, kp, nHat};
	res = pcubature(fdim, intFx, &fdata, dim, &xmin, &xmax, maxEval, 
		absError, relError, ERROR_INDIVIDUAL, &val, &err);
	std::cout << "return value: " << res << std::endl;
	std::cout << "integral value: " << val << std::endl;
	std::cout << "error value: " << err << std::endl;
	return val;
}
