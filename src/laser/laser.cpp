#include "laser.h"
#include "maths/constants.h"

#include <iostream>

Laser::Laser() {}
Laser::Laser(double peak_I0_wcm2, double waist_um, double wavelength_nm) :
    _peak_I0_wcm2(peak_I0_wcm2), _waist_um(waist_um), _wavelength_nm(wavelength_nm) {
    _freq = LnmToAU / _wavelength_nm;
    _wavelength = nmToAU * _wavelength_nm;
    _peak_I0 = Wcm2ToAU * _peak_I0_wcm2;
    _waist = umToAU * _waist_um;
    _rayleigh = PI*_waist*_waist / _wavelength;
    _porras = 0;                            // Gouy phase only
}
Laser::Laser(double peak_I0_wcm2, double waist_um, double wavelength_nm, double porras) :
    _peak_I0_wcm2(peak_I0_wcm2), _waist_um(waist_um), _wavelength_nm(wavelength_nm), _porras(porras) {
    _freq = LnmToAU / _wavelength_nm;
    _wavelength = nmToAU * _wavelength_nm;
    _peak_I0 = Wcm2ToAU * _peak_I0_wcm2;
    _waist = umToAU * _waist_um;
    _rayleigh = PI*_waist*_waist / _wavelength;
}
double Laser::IntensityAt(point3 pos) const {
    double r = sqrt(pos.x*pos.x + pos.y*pos.y);
    double z = pos.z;
    double wz = Radius(z);
    double I = _peak_I0 * (_waist/wz)*(_waist/wz) * exp(-2.*(r/wz)*(r/wz));
    // std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    // std::cout << I << std::endl;
    // exit(0);
    return I;
}
double Laser::GouyPhase(point3 pos) const {
    double r = sqrt(pos.x*pos.x + pos.y*pos.y);
    double z = pos.z;
    return -atan2(z, _rayleigh);
}
double Laser::FocalPhase(point3 pos) const {
    // g0 * [1 + (r/w(z))^2] / (z/z0 + z0/z)
    double r = sqrt(pos.x*pos.x + pos.y*pos.y);
    double z = pos.z;
    return _porras * (1 - 2*(r*r/Radius(z)/Radius(z))) / (Curvature(z) / _rayleigh);
}
double Laser::Phase(point3 pos) const {
    double phase = 0;
    double r = sqrt(pos.x*pos.x + pos.y*pos.y);
    double z = pos.z;
    return 0;
    if (_porras == -1)   // not so good
        return phase;
    
    phase += GouyPhase(pos);
    if (_porras == 0) 
        return phase;

    return phase + FocalPhase(pos);
}

double Laser::PeakI0() const {
    return _peak_I0;
}
double Laser::Waist() const {
    return _waist;
}
double Laser::Radius(double z) const {
    return _waist * sqrt(1 + (z/_rayleigh) * (z/_rayleigh));
}
double Laser::Freq() const {
    return _freq;
}
double Laser::Wavelength() const {
    return _wavelength;
}
double Laser::RayleighLength() const {
    return _rayleigh;
}
double Laser::Porras() const {
    return _porras;
}
double Laser::Curvature(double z) const {
    return _rayleigh * (z / _rayleigh + _rayleigh / z);
}