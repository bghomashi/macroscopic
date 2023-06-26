#include "gas_jet.h"
#include "maths/constants.h"


CylindricalGasJet::CylindricalGasJet() : gen(rd()), dis_xz(0., 1.), dis_y(-_length/2., _length/2.), dis_t(0., 2.*PI), _density(0), _sigma(0) {}
CylindricalGasJet::CylindricalGasJet(double density, double sigma, double radius, double length, size_t n) : 
_n(n), _length(length), _radius(radius), _sigma(sigma), _density(density) {
}



void CylindricalGasJet::SampleCylinder(Laser& laser) {
    if (_n == 1) {
        double R = 0;
        double t = 0;
        double z = 0;
        double x = 0;
        double y = 0;
        double d = 1;
        double I = laser.IntensityAt(sqrt(x*x + y*y), z);
        double P = laser.Phase(sqrt(x*x + y*y), z);
        _cells.push_back({x, y, z, d, I, P});
    } else {
        for (int i = 0; i < _n;) {
            double R = _radius*sqrt(dis_xz(gen)); //gives uniform samplying of the circular cross section
            double t = dis_t(gen);
            double z = R*cos(t);
            double x = R*sin(t);
            double y = dis_y(gen);
            double d = 1;//_density*exp(-0.5*R*R/_sigma/_sigma);
            double I = laser.IntensityAt(sqrt(x*x + y*y), z);
            double P = laser.Phase(sqrt(x*x + y*y), z);

            if (I/laser.PeakI0() > 0.032) {
                _cells.push_back({x, y, z, d, I, P});
                i++;
            }
        }
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