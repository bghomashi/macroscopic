#pragma once

#include "maths/constants.h"
#include "laser/laser.h"
#include "spectrum_matrix/spectrum_matrix.h"
#include "gas_jet/gas_jet.h"

class Macroscopic {
    Laser _laser;
    SpectrumMatrix _microscopic_data;
    CylindricalGasJet _gas_jet;

    double _jetsig, _jetsig_um;
    double _density, _density_cm3;
    
    cvector MicroSpectrum(const point3& position) const;
    complex E(
            const point3& rj,                           // dipole position
            const double theta_d, const double phi_d,   // detector location
            double w, complex aj) const;                // frequency and radiation coefficient) const;

public:
    void Initialize(
        double wavelength_nm, double beam_waist, double peakI0_wcm2, 
        double gas_radius, double gas_length, size_t number_of_cells, 
        double jetsig_um, double density_cm3,
        const std::string& filename);
   
    cvector Spectrum(double td, double pd) const;
    cvector Spectrum() const;
    dvector Frequencies() const;
};
