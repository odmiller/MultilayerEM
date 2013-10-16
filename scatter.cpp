
#include <iostream>
#include <complex>
#include "materials.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);
const cdouble II = cdouble(0.0,1.0);

// From Francoeur et. al., Jour. of Quant. Spectr. & Rad. Trans. 110, 2002 (2009)
//
// Assume layers are ordered from 0 through N (left-to-right, l.e. increasing z)
// zl = z at interface between layer l-1 and layer l (l.e. left bndry of l)
//
// eps / t = permittivity / thickness of each layer
// size(eps) = N+1
// size(t) = N-1 (two outer layers extend to +-infinity)
//
// [Al; B0] = [S11 S12; S21 S22]*[A0; Bl], S = Sl
// return s = [s11-s22(0); s11-s22(1); ... s11-s22(N)]

// Fresnel coefficients (Sipe's notation)
// NOTE: kp should be normalized by k0! pass kp/k0
void reflTrans(cdouble eps0, cdouble eps1, double kp, 
	cdouble &rs, cdouble &rp, cdouble &ts, cdouble &tp) {
	cdouble w0, w1;
	w0 = sqrt( eps0 - kp*kp );
	w1 = sqrt( eps1 - kp*kp );
	rs = (w0-w1) / (w0+w1);
	rp = (w0*eps1 - w1*eps0) / (w0*eps1 + w1*eps0);
	ts = 2.*w0 / (w0+w1);
	tp = 2.*sqrt(eps0)*sqrt(eps1)*w0 / (w0*eps1 + w1*eps0);
}

// scattering matrix
void smatrix(const epsfn *eps, const double *t, int N, double w, double kp, cdouble *sTE, cdouble *sTM) {
	sTE[0] = 1;
	sTE[1] = 0;
	sTE[2] = 0;
	sTE[3] = 1;
	sTM[0] = 1;
	sTM[1] = 0;
	sTM[2] = 0;
	sTM[3] = 1;
	double tl, k0 = w/3e8;
	cdouble kl, kzl, exp1, exp2, rs, rp, ts, tp;
	int s11l, s12l, s21l, s22l;
	for(int l=1; l<N+1; ++l) {
		tl = (l>1) ? t[l-2] : 0;
		kl = sqrt(eps[l-1](w)) * k0;
		kzl = sqrt( kl*kl - kp*kp );
		exp1 = exp(II*kzl*tl);
		exp2 = exp(2.*II*kzl*tl);
		std::cout << "permittivity of layer " << l << " = " << eps[l](w) << std::endl;
		reflTrans( eps[l-1](w), eps[l](w), kp/k0, rs, rp, ts, tp);
		std::cout <<"l: " << l << "  rs: " << rs << "  rp: " << rp << "  ts: " 
			<< ts << "  tp: " << tp << "  exp1: " << exp1 << "  exp2: " << exp2 
			<< "  kl: " << kl << "  kp: " << kp << "  kzl: " << kzl << std::endl;

		s11l = 4*l; // subtract 4 to get previous iteration's entries
		s12l = 4*l+1;
		s21l = 4*l+2;
		s22l = 4*l+3;

		sTE[s11l] = ( sTE[s11l-4] * ts * exp1 )
				/ ( 1. - sTE[s12l-4] * rs * exp2 ); //s11
		sTE[s12l] = ( sTE[s12l-4] * exp2 - rs )
				/ ( 1. - sTE[s12l-4] * rs * exp2 );
		sTE[s21l] = sTE[s11l] * sTE[s22l-4] * rs * exp1 / ts + sTE[s21l-4];
		sTE[s22l] = sTE[s22l-4] * ( rs * sTE[s12l] + 1. ) * exp1 / ts;

		sTM[s11l] = ( sTM[s11l-4] * tp * exp1 )
				/ ( 1. - sTE[s12l-4] * rp * exp2 ); //s11
		sTM[s12l] = ( sTM[s12l-4] * exp2 - rp )
				/ ( 1. - sTM[s12l-4] * rp * exp2 );
		sTM[s21l] = sTM[s11l] * sTM[s22l-4] * rp * exp1 / tp + sTM[s21l-4];
		sTM[s22l] = sTM[s22l-4] * ( rp * sTM[s12l] + 1. ) * exp1 / tp;
	}
}

/*
// Si at 400nm for testing
cdouble epsSi(double w) {
	cdouble n = cdouble(5.57, 0.387);
	return n*n;
}

int main() {
	// test Fresnel coefficients
	cdouble rs, rp, ts, tp;
	cdouble eps0 = 1.0;
	cdouble eps1 = 3.5*3.5;
	reflTrans(eps0, eps1, sin(M_PI/6.), rs, rp, ts, tp);
	std::cout << "rs: " << rs << "\nrp: " << rp << "\nts: " 
		<< ts << "\ntp: " << tp << "\n" << std::endl;

	// test scattering matrix
	double w = 2*M_PI*3e8/400e-9;
	int noLayers = 5;
	epsfn eps[] = {epsVac, epsSi, epsVac, epsSi, epsVac};
	double t[] = { 1.0e-6, 0.5e-6, 0.25e-6 };
	cdouble *sTE = new cdouble[4*noLayers], *sTM = new cdouble[4*noLayers];
	smatrix(eps, t, noLayers-1, w, w/3e8 * sin(M_PI/3), sTE, sTM); 
	for(int i=0; i<4*noLayers; ++i)
		std::cout << "sTE[" << i << "] = " << sTE[i] << std::endl;
	std::cout << std::endl;
	for(int i=0; i<4*noLayers; ++i)
		std::cout << "sTM[" << i << "] = " << sTM[i] << std::endl;
}
*/

