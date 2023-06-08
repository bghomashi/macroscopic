#pragma once 

#include <complex>
#include <vector>
#include <map>
#include "maths/constants.h"
// #include "spline.h"

class SpectrumMatrix {
    dvector frequencies, intensities;
    cvector data;
    // std::vector<Spline::Cubic> splines;

    complex& get(int i, int j);
    complex get(int i, int j) const;
    complex lerp(double x0, double x1, complex y0, complex y1, double x) const;
    
public:
    cvector Lerp(double intensity) const;

    const dvector& Frequencies() const;
    static SpectrumMatrix Load(const std::string& filename);
};
