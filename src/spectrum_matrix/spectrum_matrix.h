#pragma once 

#include <complex>
#include <vector>
#include <map>
#include "maths/constants.h"
#include "spline.h"

class SpectrumMatrix {
    dvector frequencies, intensities;
    cvector data;
    std::vector<Spline::Cubic> splines;

    complex lerp(double x0, double x1, complex y0, complex y1, double x) const;
    
public:
    complex& get(int i, int f);
    complex get(int i, int f) const;
    cvector Lerp(double intensity) const;

    const dvector& Frequencies() const;
    const dvector& Intensities() const;
    static SpectrumMatrix Load(const std::string& filename);
};
