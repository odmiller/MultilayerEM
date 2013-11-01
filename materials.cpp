
#include <iostream> // for printEps
#include <complex>

typedef std::complex<double> cdouble;
typedef cdouble (*epsfn)(double);
const cdouble II(0.0,1.0);

double c0 = 299792458;
double q = 1.602176565e-19;
double eps0 = 8.85418782e-12;
double mu0 = M_PI*4e-7;
double h = 6.62606957e-34;
double hbar = h / (2.*M_PI);
double wToEV = hbar / q;

// add logspace option
void printEps(epsfn fx, double w1, double w2, int Nw) {
	double dw = (Nw>1) ? (w2 - w1) / (Nw - 1) : 2*(w2-w1);
	for(double w=w1; w<=w2; w += dw) 
		std::cout << w << " " << fx(w) << std::endl;
}

cdouble epsAu(double w) {
	double epsInf = 1., wp = 1.3594122e16, g = 1.0495462e14;
	return epsInf - wp*wp / (w*w + II*g*w);
}

// for now, just a simple Drude model
// could also have a multi-oscillator model
cdouble epsAg(double w) {
	double epsInf = 1., wp = 1.3689e16, g = 2.7347e13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

// ref: Joulain et. al. PRB 68, 245405 (2003)
cdouble epsAl(double w) {
	double epsInf = 1., wp = 1.747e16, g = 7.596e13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

cdouble epsCr(double w) {
	double epsInf = 1., wp = 7.1942e14, g = 7.6921e13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

cdouble epsW(double w) {
	double epsInf = 1., wp = 8.8118e15, g = 7.5963e13;
	return epsInf - wp*wp / (w*w + II*g*w);
}

// SiC, from Ben-Abdallah et. al. J. Appl. Phys. 106, 044306 (2009)
cdouble epsSiC(double w) {
	double epsInf, wLO, wTO, g;
	epsInf = 6.7;
	wLO = 18.253e13;
	wTO = 14.937e13;
	g = 8.966e11;
	return epsInf * (1. + (wLO*wLO - wTO*wTO) / (wTO*wTO - w*w - II*g*w)); // same as in epsCBN
}

cdouble epsSiDoped(double w, double N) {
	double epsf = 1.0035;
	double eps0Si = 11.87;
	double w0 = 6.6e15;
	double epsSi = epsf + (eps0Si - epsf) * w0*w0 / (w0*w0 - w*w);
	if(N==0)
		return epsSi;
	N = N * 1e6; // convert to /m^3
	double q = 1.602176565e-19;
	double rho = 1.1e-4;
	double me = 9.10938215e-31;
	double m = 0.34*me;
	double eps0 = 8.85418782e-12;
	double wp = sqrt(N * q * q/(eps0 * m));
	double g = N * q * q * rho / m;
	return epsSi - wp*wp / (w*w + II*g*w);
}

// c-BN, from Francoeur et. al. JQSRT 110, 2002 (2009)
cdouble epsCBN(double w) {
	double epsInf, wLO, wTO, g;
	epsInf = 4.46;
	wLO = 2.451e14;
	wTO = 1.985e14;
	g = 9.934e11;
	return epsInf * (w*w - wLO*wLO + II*g*w) / (w*w - wTO*wTO + II*g*w);
}

cdouble epsVac(double w) {
	return cdouble(1.,0.);
}

// for next three methods:
//   with c++11, can use lambda expression to change signature
//   e.g. auto f = [=] (double d) { return epsCont(d,1); };
cdouble epsConst(double w, cdouble eps) {
	return eps;
}

cdouble epsDrude(double w, double wp, double g) {
	return cdouble( 1 - wp*wp /( w*w + g*g ), g*wp*wp /( w*w*w + g*g*w ) );
}

cdouble epsDrude(double w, double wp, double g, double epsInf) {
	return cdouble( epsInf - wp*wp /( w*w + g*g ), g*wp*wp /( w*w*w + g*g*w ) );
}
