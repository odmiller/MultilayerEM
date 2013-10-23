
/* This is exactly the same as SMatrix.cpp, except with a slower recursion 
 * formula used to compute the SMatrices.  Used to verify formulas in 
 * SMatrix.cpp
 */

#include <iostream>
#include <complex>
#include "SMatrix.hpp"
#include "SMatrixF.hpp"

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);
const cdouble II(0.0,1.0);

// fwd decl's of some utility fxn's
void eyeS(cdouble *S, int i);
int ind11(int a, int b, int N);
int ind12(int a, int b, int N);
int ind21(int a, int b, int N);
int ind22(int a, int b, int N);

SMatrixF::SMatrixF(const mlgeo *g, double k0In, double kpIn)
	: kz(new cdouble[g->N+1]), k(new cdouble[g->N+1]), kp(kpIn), k0(k0In), N(g->N),
	sTE(new cdouble[4*(N+2)*(N+1)/2]), sTM(new cdouble[4*(N+2)*(N+1)/2]) 
{
	initMatrix(g->eps, g->d);
}

// units of k0, kp, only matter relative to d (only m*d products ever used)
SMatrixF::SMatrixF(const cdouble *eps, const double *d, int NIn, double k0In, double kpIn)
	: kz(new cdouble[N+1]), k(new cdouble[N+1]), kp(kpIn), k0(k0In), N(NIn),
	sTE(new cdouble[4*(N+2)*(N+1)/2]), sTM(new cdouble[4*(N+2)*(N+1)/2]) // member init. list
{
	initMatrix(eps, d);
}

// simple destructor
SMatrixF::~SMatrixF() {
	delete[] kz;
	delete[] k;
	delete[] sTE;
	delete[] sTM;
}

void SMatrixF::initMatrix(const cdouble *eps, const double *d) {
	double dl;
	cdouble exp1, exp2, r, t, *s;
	int m, n; // going to set S(0,n) and S(m,N), m=0:N,n=0:N

	for(int pol=0; pol<2; ++pol) {
		s = (pol==TE) ? sTE : sTM;

		for(m=0; m<=N; ++m) {
			k[m] = sqrt(eps[m]) * k0;
			kz[m] = sqrt(k[m] * k[m] - kp * kp);
		}

		// S(m,n) for m=0:N, n=m:N
		for(m=0; m<=N; ++m) {
			for(n=m; n<=N; ++n) { 
				//std::cout << m << " " << n << std::endl;
				//std::cout << "ind11: " << ind11(m,n,N) << std::endl;
				if(m==n) {
					eyeS(s,ind11(m,n,N));
					continue;
				}
				dl = (n>1) ? d[n-2] : 0;
				exp1 = exp(II*kz[n-1]*dl);
				exp2 = exp(2.*II*kz[n-1]*dl);
				reflTrans( eps[n-1], eps[n], kp/k0, pol, r, t );
				
				//std::cout << "[SMatrix]: " << m << " " << n << " " << r << " " << t << std::endl;
				if(m==0) {
					std::cout << "[SMatrixF]: " << m << " " << n << " " << r << " " << t 
						<< exp1 << " " << exp2 << " " << s[ind11(m,n-1,N)] << " " 
						<< s[ind12(m,n-1,N)] << std::endl;
				}
				s[ind11(m,n,N)] = ( s[ind11(m,n-1,N)] * t * exp1 )
						/ ( 1. - s[ind12(m,n-1,N)] * r * exp2 ); 
				s[ind12(m,n,N)] = ( s[ind12(m,n-1,N)] * exp2 - r )
						/ ( 1. - s[ind12(m,n-1,N)] * r * exp2 );
				s[ind21(m,n,N)] = s[ind11(m,n,N)] * s[ind22(m,n-1,N)] * r * exp1 / t + s[ind21(m,n-1,N)];
				s[ind22(m,n,N)] = s[ind22(m,n-1,N)] * ( r * s[ind12(m,n,N)] + 1. ) * exp1 / t;
			}
		}
	}
}

cdouble SMatrixF::S11(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind11(a,b,N)];
	else
		return sTE[ind11(a,b,N)];
}

cdouble SMatrixF::S12(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind12(a,b,N)];
	else
		return sTE[ind12(a,b,N)];
}

cdouble SMatrixF::S21(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind21(a,b,N)];
	else
		return sTE[ind21(a,b,N)];
}

cdouble SMatrixF::S22(int a, int b, int pol) const {
	if(pol==TM)
		return sTM[ind22(a,b,N)];
	else
		return sTE[ind22(a,b,N)];
}

void SMatrixF::printSMatrix() {
	for(int pol=0; pol<=1; ++pol) {
		std::cout << "-------" << std::endl;
		if(pol==TE)
			std::cout << "pol: " << "TE" << std::endl;
		else
			std::cout << "pol: " << "TM" << std::endl;
		std::cout << "-------" << std::endl;
		for(int m=0; m<=N; ++m) {
			for(int n=m; n<=N; ++n) {
				std::cout << "  s11(" << m << "," << n << "): " 
					<< S11(m,n,pol) << std::endl;
				std::cout << "  s12(" << m << "," << n << "): " 
					<< S12(m,n,pol) << std::endl;
				std::cout << "  s21(" << m << "," << n << "): " 
					<< S21(m,n,pol) << std::endl;
				std::cout << "  s22(" << m << "," << n << "): " 
					<< S22(m,n,pol) << std::endl;
			}
		}
	}
}

/*
int main() {
	// test Fresnel coefficients
	cdouble rs, rp, ts, tp;
	cdouble *epsA;
	cdouble eps0 = 1.0;
	cdouble eps1 = 3.5*3.5;
	reflTrans(eps0, eps1, sin(M_PI/6.), 0, rs, ts);
	reflTrans(eps0, eps1, sin(M_PI/6.), 1, rp, tp);
	std::cout << "rs: " << rs << "\nrp: " << rp << "\nts: " 
		<< ts << "\ntp: " << tp << "\n" << std::endl;

	// test scattering matrix
	double w = 2*M_PI*3e8/400e-9;
	int noLayers = 5;
	epsfn epsSi = [=] (double d) { return epsConst(d,cdouble(30.875,4.311)); };
	epsfn eps[] = {epsVac, epsSi, epsVac, epsSi, epsVac};
	epsA = new cdouble[noLayers];
	for(int i=0; i<noLayers; ++i)
		epsA[i] = eps[i](w);
	double d[] = { 0.25e-6, 0.5e-6, 0.25e-6 };
	mlgeo *g = new mlgeo;
	initGeo(epsA, d, noLayers-1, g);
	printGeo(g);	
	SMatrix *S = new SMatrix(g, w/3.e8, w/3.e8*sqrt(21./25.));
	//SMatrix *S = new SMatrix(epsA, d, noLayers-1,w/3.e8,w/3.e8*sqrt(21./25.));
	S->printSMatrix();
}
*/

/*void eyeS(cdouble *S, int i) { // insert identity matrix
	S[i+0] = 1; S[i+1] = 0;
	S[i+2] = 0; S[i+3] = 1;
}*/
int ind11(int a, int b, int N) {
	return 4 * (a*(N+1)-a*(a-1)/2+(b-a));
}
int ind12(int a, int b, int N) {
	return 4 * (a*(N+1)-a*(a-1)/2+(b-a)) + 1;
}
int ind21(int a, int b, int N) {
	return 4 * (a*(N+1)-a*(a-1)/2+(b-a)) + 2;
}
int ind22(int a, int b, int N) {
	return 4 * (a*(N+1)-a*(a-1)/2+(b-a)) + 3;
}
