#pragma once

#include <complex>
#include <cassert>
#include <vector>

using namespace std::complex_literals;          // enables 'i' in complex literals: complex x = 2.5i;

typedef std::complex<double> complex;
typedef std::vector<complex> cvector;
typedef std::vector<double> dvector;

const double PI =       3.1415926535897932384626433832795;
const double C =        137.0;                  // speed of light in AU
const double LnmToAU =  45.5513500662;          // wavelength in nm to AU (wavelength in *denominator*)
const double nmToAU =   1./ 0.052917721092;     // length in nm to AU
const double cmToAU =   1./ 0.52917721092e-8;    // length in cm to AU
const double Wcm2ToAU = 1./3.51e16;             // 
const double umToAU =   1.89393939394e4;        // length in um to AU

struct point2 {
    double x, y, z;
};
struct point3 {
    double x, y, z;
};
