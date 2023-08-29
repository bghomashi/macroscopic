#pragma once

#include "maths/constants.h"
#include "laser/laser.h"
#include "spectrum_matrix/spectrum_matrix.h"
#include "points/points.h"

class Macroscopic {
    Laser _laser;
    SpectrumMatrix _microscopic_data;
    points _points;
    cvector MicroSpectrum(const point3& position) const;
    complex E(
            const point3& rj,                           // dipole position
            const double theta_d, const double phi_d,   // detector location
            double w, complex aj) const;                // frequency and radiation coefficient) const;

public:
    void Initialize(Laser& laser, points samples, const std::string& filename);
   
    cvector Spectrum(double td, double pd) const;
    cvector Spectrum() const;
    dvector Frequencies() const;
};
