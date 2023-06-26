#pragma once

class Laser {
    double _peak_I0, _waist, _wavelength, _freq, _rayleigh;     // in atomic units
    double _peak_I0_wcm2, _wavelength_nm, _waist_um;            // in natural units
    double _porras;                                             // unitless
public:
    Laser();
    Laser(double peak_I0_wcm2, double waist_um, double wavelength_nm);
    Laser(double peak_I0_wcm2, double waist_um, double wavelength_nm, double porras);
    
    double IntensityAt(double r, double z) const;
    double GouyPhase(double r, double z) const;
    double FocalPhase(double r, double z) const;
    double Phase(double r, double z) const;

    double PeakI0() const;                  // in atomic units
    double Waist() const;                   // in atomic units
    double Radius(double z) const;          // in atomic units
    double Curvature(double z) const;       // in atomic units
    double Freq() const;                    // in atomic units
    double Wavelength() const;              // in atomic units
    double RayleighLength() const;          // in atomic units
    double Porras() const;
};