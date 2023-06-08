#include "gas_jet.h"


CylindricalGasJet::CylindricalGasJet() : gen(rd()), dis_xz(0., 1.), dis_y(-_length/2., _length/2.), dis_t(0., 2.*M_PI), _density(0), _sigma(0) {}
CylindricalGasJet::CylindricalGasJet(double density, double sigma, double radius, double length, size_t n) : 
_n(n), _length(length), _radius(radius), _sigma(sigma), _density(density) {
    SampleCylinder();
}



void CylindricalGasJet::SampleCylinder() {
    for (int i = 0; i < _n; i++) {
        double R = _radius*sqrt(dis_xz(gen)); //gives uniform samplying of the circular cross section ??
        double t = dis_t(gen);
        double z = R*cos(t);
        double x = R*sin(t);
        double y = dis_y(gen);
        double d = _density*exp(-0.5*R*R/_sigma/_sigma);
        _cells.push_back({x, y, z, d, 0});
    }
}

CylindricalGasJet& CylindricalGasJet::operator=(const CylindricalGasJet& o) {
    _density = o._density;
    _sigma = o._sigma;
    _radius = o._radius;
    _length = o._length;
    _n = o._n;
    dis_xz = o.dis_xz;
    dis_y = o.dis_y;
    dis_t = o.dis_t;
    _cells = o._cells;

    return *this;
}