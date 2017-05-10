#ifndef _basic_types_h
#define _basic_types_h


#include <limits>

namespace nr
{
	
	// basic type names (redefine if your bit lengths don't match)

	typedef int Int; // 32 bit integer
	typedef unsigned int Uint;

	#ifdef _MSC_VER
	typedef __int64 Llong; // 64 bit integer
	typedef unsigned __int64 Ullong;
	#else
	typedef long long int Llong; // 64 bit integer
	typedef unsigned long long int Ullong;
	#endif

	typedef char Char; // 8 bit integer
	typedef unsigned char Uchar;

	typedef double Doub; // default floating type
	typedef long double Ldoub;

	typedef std::complex<double> Complex; // default complex type

	typedef bool Bool;

	// NaN: uncomment one of the following 3 methods of defining a global NaN
	// you can test by verifying that (NaN != NaN) is true

	static const Doub NaN = std::numeric_limits<Doub>::quiet_NaN();

	//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
	//double NaN = *( double* )proto_nan;

	//Doub NaN = sqrt(-1.);

}

#endif	//_basic_types_h