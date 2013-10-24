
#include <iostream>
#include <complex>
#include "cubature.h"
#include "multilayer.hpp"
#include "SMatrix.hpp"

struct geoData {
	const mlgeo *g;
	int l, s;
	double zl, k0, kp, nHat;
};

int intFxKp(unsigned ndim, const double *kp, void *fdata, unsigned fdim, double *fval) {
	geoData *f = (geoData *) fdata;
	*fval = flux(f->g, f->l, f->s, f->zl, f->k0, (*kp) * f->k0, f->nHat);
	//std::cout << "kp: " << *kp << std::endl;
	//std::cout << "fval: " << *fval << std::endl;
	return 0;
}

int intFxZ(unsigned ndim, const double *zs, void *fdata, unsigned fdim, double *fval) {
	geoData *f = (geoData *) fdata;
	*fval = flux(f->g, f->l, f->s, f->zl, f->k0, f->kp, f->nHat, *zs);
	return 0;
}

// integrate over kp (integration over emitter layer assumed)
double fluxKpInt(const mlgeo *g, int l, int s, double zl, double k0, double nHat) {
	int res;
	unsigned fdim = 1, dim = 1;
	double xmin = 0, xmax = 10000., val, err; // transform later
	double absError = 0, relError = 1e-4;
	size_t maxEval = 1e5;
	geoData fdata = {g, l, s, zl, k0, 0, nHat};
	res = hcubature(fdim, intFxKp, &fdata, dim, &xmin, &xmax, maxEval,
					absError, relError, ERROR_INDIVIDUAL, &val, &err);
	if(res!=0)
		std::cout << "k0: " << k0 << "  return value: " << res << std::endl;
	//std::cout << "return value: " << res << std::endl;
	//std::cout << "integral value: " << val << std::endl;
	//std::cout << "error value: " << err << std::endl;
	return k0*val; // in units of k0 above
}

// integrate over emitter layer (have analytical formula - should be unnecessary!)
double fluxZInt(mlgeo *g, int l, int s, double zl, double k0, double kp, double nHat) {
	int res;
	unsigned fdim = 1, dim = 1;
	double xmin = 0, xmax = g->d(s), val, err; // s cannot equal 0 or N 
	double absError = 0, relError = 1e-6;
	size_t maxEval = 1e5;
	geoData fdata = {g, l, s, zl, k0, kp, nHat};
	res = pcubature(fdim, intFxZ, &fdata, dim, &xmin, &xmax, maxEval, 
		absError, relError, ERROR_INDIVIDUAL, &val, &err);
	std::cout << "return value: " << res << std::endl;
	std::cout << "integral value: " << val << std::endl;
	std::cout << "error value: " << err << std::endl;
	return val;
}
