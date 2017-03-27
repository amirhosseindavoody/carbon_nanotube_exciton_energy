#ifndef constants_h
#define constants_h

#include <complex>
#include <math.h> // using pow

using namespace std;

namespace constants
{
	// constants();

	static const double pi=3.141592;
	
	static const double eV=pow(1.6,-19.0); //[Joules]
	static const double hb=pow(6.5,-16.0)*eV; //[eV.s]
	static const double kb=pow(1.3865,-23.0); //[Joules/Kelvin]

	static const double eps0 = pow(8.85,-12); // permittivity of free space
	static const double q0 = pow(1.6,-19); // charge of electron
};

#endif // constants_h