
#include <iostream>
#include <complex>
#include "cubature.h"
#include "fieldFunctions.hpp"
#include "SMatrix.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);
const double c0 = 299798452;

struct fluxData {
	const mlgeo &g;
	const int *s;
	const int *l;
	const double *zl;
	const double *nHat;
	int Ns, Nl;
	double k0;
};

struct fluxData2 {
	const epsfn *eps;
	const double *d;
	const int *s;
	const int *l;
	const double *zl;
	const double *nHat;
	double lscale, temp; // temp = temperature
	int numLayers, Ns, Nl;
};

struct dosData {
	const mlgeo &g;
	const int *l;
	const double *zl;
	double k0;
};

int fluxFx(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	double t0 = x[0];
	double kp = t0 / (1. - t0);
	fluxData *f = (fluxData *) fdata;
	SMatrix S = SMatrix(f->g, f->k0, kp);
	//std::cout << f->k0 << " " << kp << std::endl;
	fval[0] = 0;
	for (int il = 0; il < f->Nl; ++il)
		fval[0] += flux(f->g, S, f->l[il], f->zl[il], f->s, f->Ns, f->nHat[il]);
	if(fval[0] != fval[0]) // isnan
		fval[0] = 0;
	else
		fval[0] /= pow(1. - t0, 2);
	return 0;
}

// x = [k0; kp]
int fluxFx2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	double t0 = x[0];
	double k0 = t0 / (1. - t0);
	double t1 = x[1];
	double kp = t1 / (1. - t1);

	fluxData2 *f = (fluxData2 *) fdata;
	double w = k0 * c0 / f->lscale;
	cdouble *epsV = new cdouble[f->numLayers];
	for (int i = 0; i < f->numLayers; ++i)
		epsV[i] = f->eps[i](w);
	mlgeo g = mlgeo(epsV, f->d, f->numLayers - 1);
	SMatrix S = SMatrix(g, k0, kp);
	fval[0] = 0;
	for (int il = 0; il < f->Nl; ++il)
		fval[0] += flux(g, S, f->l[il], f->zl[il], f->s, f->Ns, f->nHat[il]);
	if( fval[0] != fval[0] ) // isnan
		fval[0] = 0.;
	else
		fval[0] *= mean_energy(k0, f->lscale, f->temp) / pow((1. - t0) * (1. - t1), 2);
	//std::cout << t0 << " " << t1 << " " << k0 << " " << kp/k0 << " " << *fval 
			//<< " " << mean_energy(k0, f->lscale, f->temp) << std::endl;
	return 0;
}

int dosFx(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	double t0 = x[0];
	double kp = t0 / (1. - t0);
	dosData *f = (dosData *) fdata;
	SMatrix S = SMatrix(f->g, f->k0, kp);
	for (unsigned il = 0; il < fdim; ++il) {
		fval[il] = dos(S, f->l[il], f->zl[il]) / pow(1. - t0, 2);
		if(fval[il] != fval[il])
			fval[il] = 0.;
	}
	//std::cout << "k0: " << f->k0 << "  kp: " << kp << " " << *fval << std::endl;
	return 0;
}

// integrate over kp and w (integration over emitter layer assumed)
// This allows Nl flux planes (each with a corresponding zl and nHat) given by l
double fluxKpWInt(const epsfn *eps, const double *d, const int *l, const double *zl, 
				   const double *nHat, const int *s, int numLayers, int Nl, int Ns, 
				   double lscale, double temp) {
	unsigned fdim = 1, dim = 2;
	double xmin[2] = {0,0};
	double xmax[2] = {1,1.};
	double val, err;
	double absError = 0, relError = 1e-6;
	size_t maxEval = 1e5;
	fluxData2 fdata = {eps, d, s, l, zl, nHat, lscale, temp, numLayers, Ns, Nl};
	int res = hcubature(fdim, fluxFx2, &fdata, dim, xmin, xmax, maxEval,
					absError, relError, ERROR_INDIVIDUAL, &val, &err);
	if (res != 0)
		std::cout << "Return value: " << res << std::endl;
	return val; // Note that the integral is over normalized kp/k0
}

// integrate over kp (integration over emitter layer assumed)
// This allows Nl flux planes (each with a corresponding zl and nHat) given by l
double fluxKpInt(const mlgeo &g, const int *l, const double *zl, const double *nHat, 
				 const int *s, int Nl, int Ns, double k0) {
	unsigned fdim = 1, dim = 1;
	double xmin = 0, xmax = 1, val, err; // transform later
	double absError = 0, relError = 1e-4;
	size_t maxEval = 1e5;
	fluxData fdata = {g, s, l, zl, nHat, Ns, Nl, k0};
	int res = hcubature(fdim, fluxFx, &fdata, dim, &xmin, &xmax, maxEval,
						absError, relError, ERROR_INDIVIDUAL, &val, &err);
	if (res != 0)
		std::cout << "k0: " << k0 << "  return value: " << res << std::endl;
	return val;
}

double *dosKpInt(const mlgeo &g, const int *l, const double *zl, int Nl, double k0) {
	int res;
	unsigned fdim = Nl, dim = 1;
	double xmin = 0, xmax = 1;
	double *val = new double[fdim];
	double *err = new double[fdim]; // transform later
	double absError = 0, relError = 1e-4;
	size_t maxEval = 1e5;
	dosData fdata = {g, l, zl, k0};
	res = hcubature(fdim, dosFx, &fdata, dim, &xmin, &xmax, maxEval,
					absError, relError, ERROR_INDIVIDUAL, val, err);
	if (res != 0)
		std::cout << "k0: " << k0 << "  return value: " << res << std::endl;
	delete[] err;
	return val;
}

