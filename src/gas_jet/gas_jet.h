#pragma once

#include <random>
#include "maths/constants.h"
#include "laser/laser.h"

class CylindricalGasJet {
    double _density, _sigma;
    double _radius, _length;
    size_t _n;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dis_xz;
    std::uniform_real_distribution<double> dis_y;
    std::uniform_real_distribution<double> dis_t;
public:
    struct Cell {
        point3 pos;
        double density;
        double intensity;
        double phase;
    };
    std::vector<Cell> _cells;

    CylindricalGasJet();
    CylindricalGasJet(double density, double sigma, double radius, double length, size_t num);

    void SampleCylinder(Laser& laser);
    CylindricalGasJet& operator=(const CylindricalGasJet& o);
};

