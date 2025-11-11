#ifndef CHAPLIN_H
#define CHAPLIN_H

#include <iostream>

/// @todo Decide on how to handle complex numbers
// Fortran's COMPLEX*16 layout:
struct mycomplex {
    double re;
    double im;
};

/// functions hpli_ are defined in libchaplin.a
extern "C"	{
	mycomplex hpl1_(int* n1, mycomplex* x);
	mycomplex hpl2_(int* n1, int* n2, mycomplex* x);
	mycomplex hpl3_(int* n1, int* n2, int* n3, mycomplex* x);
	mycomplex hpl4_(int* n1, int* n2, int* n3, int* n4, mycomplex* x);
}

/// @brief generic interface for C++ to the Fortran routine
inline double HPL1(int n1, double x) { 
	mycomplex xc{x,0.0};
	double result(0.0);
	/// Casting everything to real numbers might be dangerous!!!
	/// I am doing it because I expect all imaginary parts to cancel in the computation
	/// of coefficient functions and there are no products of HPLs, so there are no
	/// additional real valued contributions from products of imaginary parts.
	result += hpl1_(&n1, &xc).re;
	// if(result == NAN) std::cout << "NaN in HPL1("<<n1<<","<<x<<")";
	return result;
};

/// @brief generic interface for C++ to the Fortran routine
inline double HPL2(int n1, int n2, double x) { 
	mycomplex xc{x,0.0};
	double result(0.0);
	/// Casting everything to real numbers might be dangerous!!!
	/// I am doing it because I expect all imaginary parts to cancel in the computation
	/// of coefficient functions and there are no products of HPLs, so there are no
	/// additional real valued contributions from products of imaginary parts.
	result += hpl2_(&n1, &n2, &xc).re;
	// if(result == NAN) std::cout << "NaN in HPL2("<<n1<<","<<n2<<","<<x<<")";
	return result;
};

/// @brief generic interface for C++ to the Fortran routine
inline double HPL3(int n1, int n2, int n3, double x) { 
	mycomplex xc{x,0.0};
	double result(0.0);
	/// Casting everything to real numbers might be dangerous!!!
	/// I am doing it because I expect all imaginary parts to cancel in the computation
	/// of coefficient functions and there are no products of HPLs, so there are no
	/// additional real valued contributions from products of imaginary parts.
	result += hpl3_(&n1, &n2, &n3, &xc).re;
	// if(result == NAN) std::cout << "NaN in HPL3("<<n1<<","<<n2<<","<<n3<<","<<x<<")";
	return result;
};

/// @brief generic interface for C++ to the Fortran routine
inline double HPL4(int n1, int n2, int n3, int n4, double x) { 
	mycomplex xc{x,0.0};
	double result(0.0);
	/// Casting everything to real numbers might be dangerous!!!
	/// I am doing it because I expect all imaginary parts to cancel in the computation
	/// of coefficient functions and there are no products of HPLs, so there are no
	/// additional real valued contributions from products of imaginary parts.
	result += hpl4_(&n1, &n2, &n3, &n4, &xc).re;
	// if(result == NAN) std::cout << "NaN in HPL4("<<n1<<","<<n2<<","<<n3<<","<<n4<<","<<x<<")";
	return result;
};

#endif