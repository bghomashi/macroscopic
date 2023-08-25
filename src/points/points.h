#include "maths/constants.h"
#include "laser/laser.h"

class points {
public:
    struct Cell {
        point3 pos;
        double density;
        double intensity;
        double phase;
    };
    std::vector<Cell> _cells;

    double cutoff;


    points();
    points(std::vector<point3> positions); //if you only know the positions of the points
    points(std::vector<point3> positions, Laser& laser, double cut); // if you know the positions of the points and the laser object
    void CylindricalGasJet(double density, double sigma, double radius, double length, size_t num);
    void CalcCells(Laser& laser);
};