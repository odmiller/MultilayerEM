
#ifndef MATERIALS_H
#define MATERIALS_H

#include <complex>

void printEps(std::complex<double> (*fx)(double), double w1, double w2, int Nw);
std::complex<double> epsAu(double w);
std::complex<double> epsAg(double w);
std::complex<double> epsCr(double w);
std::complex<double> epsW(double w);
std::complex<double> epsSiC(double w); 
std::complex<double> epsSiDoped(double w, double N=0.); 
std::complex<double> epsCBN(double w); // c-BN
std::complex<double> epsVac(double w);
std::complex<double> epsConst(double w, std::complex<double> eps);
std::complex<double> epsDrude(double w, double wp, double gamma); // epsInf=1
std::complex<double> epsDrude(double w, double wp, double gamma, double epsInf);

#endif
