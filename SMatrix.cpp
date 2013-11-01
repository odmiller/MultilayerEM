
#include <iostream>
#include <complex>
#include "mlgeo.hpp"
#include "SMatrix.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);
const cdouble II(0.0,1.0);

// Fresnel coefficients (Sipe's notation)
// NOTE: kp should be normalized by k0! pass kp/k0
void reflTrans(cdouble eps0, cdouble eps1, double kp, int pol,
	cdouble &r, cdouble &t) {
	cdouble w0, w1;
	w0 = sqrt( eps0 - kp*kp );
	w1 = sqrt( eps1 - kp*kp );
	if(pol==TM) {
		r = (w0*eps1 - w1*eps0) / (w0*eps1 + w1*eps0);
		t = 2.*sqrt(eps0)*sqrt(eps1)*w0 / (w0*eps1 + w1*eps0);
	} else {
		r = (w0-w1) / (w0+w1);
		t = 2.*w0 / (w0+w1);
	} 
}

// fwd decl's of some utility fxn's
void eyeS(cdouble *S, int i);
int ind11(int a, int b);
int ind12(int a, int b);
int ind21(int a, int b);
int ind22(int a, int b);

SMatrix::SMatrix(const mlgeo &g, double k0In, double kpIn)
	: kz(new cdouble[g.N+1]), k(new cdouble[g.N+1]), kp(kpIn), k0(k0In), N(g.N),
	sTE(new cdouble[4*(2*g.N+1)]), sTM(new cdouble[4*(2*g.N+1)]) // member init. list
{
	initMatrix(g);
}

// units of k0, kp, only matter relative to d (only m*d products ever used)
SMatrix::SMatrix(const cdouble *eps, const double *d, int NIn, double k0In, double kpIn)
	: kz(new cdouble[N+1]), k(new cdouble[N+1]), kp(kpIn), k0(k0In), N(NIn),
	sTE(new cdouble[4*(2*N+1)]), sTM(new cdouble[4*(2*N+1)]) // member init. list
{
	mlgeo g = mlgeo(eps, d, N);
	initMatrix(g);
}

// simple destructor
SMatrix::~SMatrix() {
	delete[] kz;
	delete[] k;
	delete[] sTE;
	delete[] sTM;
}

void SMatrix::initMatrix(const mlgeo &g) {
	double dl;
	cdouble exp1, exp2, r, t, *s;
	int m, n; // going to set S(0,n) and S(m,N), m=0:N,n=0:N

	for(int pol=0; pol<2; ++pol) {
		s = (pol==TE) ? sTE : sTM;

		// S(0,n) for n=0:N
		m = 0;
		n = 0;
		eyeS(s,ind11(m,n));
		k[0] = sqrt( g.eps(0) ) * k0;
		kz[0] = sqrt( k[0]*k[0] - kp * kp );
		for(n=1; n<=N; ++n) { 
			dl = (n>1) ? g.d(n-1) : 0;
			k[n] = sqrt( g.eps(n) ) * k0;
			kz[n] = sqrt(k[n] * k[n] - kp * kp);
			exp1 = exp(II * kz[n-1] * dl);
			exp2 = exp(2. * II * kz[n-1] * dl);
			reflTrans( g.eps(n-1), g.eps(n), kp/k0, pol, r, t );
			
			//std::cout << "[SMatrix]: " << m << " " << n << " " << r << " " << t 
				//<< exp1 << " " << exp2 << " " << s[ind11(m,n-1)] << " " 
				//<< s[ind12(m,n-1)] << std::endl;
			s[ind11(m,n)] = ( s[ind11(m,n-1)] * t * exp1 )
					/ ( 1. - s[ind12(m,n-1)] * r * exp2 ); 
			s[ind12(m,n)] = ( s[ind12(m,n-1)] * exp2 - r )
					/ ( 1. - s[ind12(m,n-1)] * r * exp2 );
			s[ind21(m,n)] = s[ind11(m,n)] * s[ind22(m,n-1)] * r * exp1 / t + s[ind21(m,n-1)];
			s[ind22(m,n)] = s[ind22(m,n-1)] * ( r * s[ind12(m,n)] + 1. ) * exp1 / t;
		}
		
		// S(m,N) for m=N:0 (backwards recurrence relations)
		m = N;
		n = N;
		eyeS(s,ind11(m,n));
		for(m=N; m>1; --m) {  // m-1 is the new layer, from m
			dl = (m>1) ? g.d(m-1) : 0;
			exp1 = exp(II*kz[m-1]*dl);
			exp2 = exp(2.*II*kz[m-1]*dl);
			reflTrans( g.eps(m), g.eps(m-1), kp/k0, pol, r, t );
			
			s[ind11(m-1,n)] = s[ind11(m,n)] * exp1 * (1. - r*r) /( t*(1.-r*s[ind21(m,n)]) );
			s[ind12(m-1,n)] = s[ind12(m,n)] + r * s[ind11(m,n)] * s[ind22(m,n)] 
							/ ( 1. - r*s[ind21(m,n)] );
			s[ind21(m-1,n)] = exp2 * (s[ind21(m,n)] - r) / (1. - r*s[ind21(m,n)]);
			s[ind22(m-1,n)] = t * exp1 * s[ind22(m,n)] / ( 1. - r*s[ind21(m,n)] );
		}
	}
}

cdouble SMatrix::S11(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind11(a,b)];
	else
		return sTE[ind11(a,b)];
}

cdouble SMatrix::S12(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind12(a,b)];
	else
		return sTE[ind12(a,b)];
}

cdouble SMatrix::S21(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind21(a,b)];
	else
		return sTE[ind21(a,b)];
}

cdouble SMatrix::S22(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind22(a,b)];
	else
		return sTE[ind22(a,b)];
}

void SMatrix::print() {
	int ind0, ind1;
	for(int pol=0; pol<=1; ++pol) {
		std::cout << "-------" << std::endl;
		if(pol==TE)
			std::cout << "pol: " << "TE" << std::endl;
		else
			std::cout << "pol: " << "TM" << std::endl;
		std::cout << "-------" << std::endl;
		for(int ind=0; ind<=1; ++ind) {
			for(int i=0; i<=N; ++i) {
				ind0 = (ind==0) ? 0 : i;
				ind1 = (ind==0) ? i : N;
				std::cout << "layer: " << i << std::endl;
				std::cout << "  s11(" << ind0 << "," << ind1 << "): " 
					<< S11(ind0,ind1,pol) << std::endl;
				std::cout << "  s12(" << ind0 << "," << ind1 << "): " 
					<< S12(ind0,ind1,pol) << std::endl;
				std::cout << "  s21(" << ind0 << "," << ind1 << "): " 
					<< S21(ind0,ind1,pol) << std::endl;
				std::cout << "  s22(" << ind0 << "," << ind1 << "): " 
					<< S22(ind0,ind1,pol) << std::endl;
			}
		}
	}
}

void eyeS(cdouble *S, int i) { // insert identity matrix
	S[i+0] = 1; S[i+1] = 0;
	S[i+2] = 0; S[i+3] = 1;
}
int ind11(int a, int b) {
	return 4*(a+b);
}
int ind12(int a, int b) {
	return 4*(a+b)+1;
}
int ind21(int a, int b) {
	return 4*(a+b)+2;
}
int ind22(int a, int b) {
	return 4*(a+b)+3;
}
