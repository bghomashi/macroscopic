#include "points.h"
#include "maths/constants.h"
#include <random>


points::points() : cutoff(0) {}

points::points(std::vector<point3> positions) : cutoff(0) {
    for (auto& pos : positions) {
        _cells.push_back({pos, 0, 0, 0});
    }
}

points::points(std::vector<point3> positions, Laser& laser, double cut) : cutoff(cut) {
    for (auto& pos : positions) {
        double I = laser.IntensityAt(pos);
        double P = laser.Phase(pos);
        if (I/laser.PeakI0() > cutoff){
            _cells.push_back({pos, 0, I, P});
        } 
    }
}

void points::CylindricalGasJet(double density, double sigma, double radius, double length, size_t num){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dis_xz(0., sigma);
    std::uniform_real_distribution<double> dis_y(-length/2., length/2.);
    if (num == 1) //if only one atom place it at the center of the jet
    {
        _cells.push_back({{0, 0, 0}, 1, 0, 0});
    }
    else
    {
        for (int i = 0; i < num;)
        {
            double x = dis_xz(gen);
            double y = dis_y(gen);
            double z = dis_xz(gen);
            double d = density*exp(-0.5*(x*x + z*z)/sigma/sigma);
            _cells.push_back({{x, y, z}, d, 0, 0});
            i++;
        }
    }
}

/*
void points::CalcCells(Laser& laser) {
    for (auto& cell : _cells) {
        cell.intensity = laser.IntensityAt(cell.pos);
        cell.phase = laser.Phase(cell.pos);
        if(cell.intensity/laser.PeakI0() < cutoff){
            _cells.erase(cell);
    }
}
*/

void points::CalcCells(Laser& laser) {
    for (int i = 0; i < _cells.size(); i++) {
        _cells[i].intensity = laser.IntensityAt(_cells[i].pos);
        _cells[i].phase = laser.Phase(_cells[i].pos);
        if(_cells[i].intensity/laser.PeakI0() < cutoff){
            _cells.erase(_cells.begin() + i);
            i--;
        }
    }
}


