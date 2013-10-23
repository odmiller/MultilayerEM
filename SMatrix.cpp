
#include <iostream>
#include <complex>
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

SMatrix::SMatrix(const mlgeo *g, double k0In, double kpIn)
	: kz(new cdouble[g->N+1]), k(new cdouble[g->N+1]), kp(kpIn), k0(k0In), N(g->N),
	sTE(new cdouble[4*(2*g->N+1)]), sTM(new cdouble[4*(2*g->N+1)]) // member init. list
{
	initMatrix(g->eps, g->d);
}

// units of k0, kp, only matter relative to d (only m*d products ever used)
SMatrix::SMatrix(const cdouble *eps, const double *d, int NIn, double k0In, double kpIn)
	: kz(new cdouble[N+1]), k(new cdouble[N+1]), kp(kpIn), k0(k0In), N(NIn),
	sTE(new cdouble[4*(2*N+1)]), sTM(new cdouble[4*(2*N+1)]) // member init. list
{
	initMatrix(eps, d);
}

// simple destructor
SMatrix::~SMatrix() {
	delete[] kz;
	delete[] k;
	delete[] sTE;
	delete[] sTM;
}

void SMatrix::initMatrix(const cdouble *eps, const double *d) {
	double dl;
	cdouble exp1, exp2, r, t, *s;
	int m, n; // going to set S(0,n) and S(m,N), m=0:N,n=0:N

	for(int pol=0; pol<2; ++pol) {
		s = (pol==TE) ? sTE : sTM;

		// S(0,n) for n=0:N
		m = 0;
		n = 0;
		eyeS(s,ind11(m,n));
		k[0] = sqrt(eps[0]) * k0;
		kz[0] = sqrt( k[0]*k[0] - kp * kp );
		for(n=1; n<=N; ++n) { 
			dl = (n>1) ? d[n-2] : 0;
			k[n] = sqrt(eps[n]) * k0;
			kz[n] = sqrt(k[n] * k[n] - kp * kp);
			exp1 = exp(II*kz[n-1]*dl);
			exp2 = exp(2.*II*kz[n-1]*dl);
			reflTrans( eps[n-1], eps[n], kp/k0, pol, r, t );
			
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
			dl = (m>1) ? d[m-2] : 0;
			exp1 = exp(II*kz[m-1]*dl);
			exp2 = exp(2.*II*kz[m-1]*dl);
			reflTrans( eps[m], eps[m-1], kp/k0, pol, r, t );
			
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

void SMatrix::printSMatrix() {
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

mlgeo::~mlgeo() {
	delete[] d;
	delete[] eps;
	delete[] z;
}

void initGeo(const cdouble *epsIn, const double *dIn, int N, mlgeo *g) {
	// create new copies
	cdouble *eps = new cdouble[N+1];
	double *d = new double[N-1];
	double *z = new double[N];
	z[0] = 0.;
	for(int i=0; i<=N; ++i) {
		eps[i] = epsIn[i];
		if(i<N-1)
			d[i] = dIn[i];
		if(i>0 && i<N)
			z[i] = z[i-1] + d[i-1];
	}
	g->eps = eps;
	g->d = d;
	g->z = z;
	g->N = N;
}

void printGeo(mlgeo *g) {
	std::cout << "-------------------------" << std::endl;
	std::cout << "Number of interfaces: " << g->N << std::endl;
	std::cout << "Number of layers: " << g->N + 1 << std::endl;
	for(int i=0; i<=g->N; ++i) {
		std::cout << "layer: " << i << std::endl;
		std::cout << "  eps: " << g->eps[i] << std::endl;
		if(i>0 && i<g->N)
			std::cout << "  d: " << g->d[i-1] << std::endl;
		if(i>0)
			std::cout << "  z: " << g->z[i-1] << std::endl;
	}
	std::cout << "-------------------------" << std::endl;
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
