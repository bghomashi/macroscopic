#pragma once

#include <complex>
#include <vector>
#include <functional>
#include <memory>
#include "maths.h"
#include "../maths/constants.h"

using namespace std::complex_literals;
typedef std::complex<double> complex;

struct Pulse {
    typedef std::shared_ptr<Pulse> Ptr_t;
    virtual point3 A(double t) const = 0;
    virtual point3 E(double t) const = 0;
    
    bool Store(const std::string& filename, const std::vector<double>& time) const;
};

struct SFA {
    double Ip;

    std::vector<Pulse::Ptr_t> pulse;
    std::vector<point3> Atot;
    std::vector<point3> integralA;
    std::vector<point3> integralAsq;
    std::vector<point3> Etot; 
    std::function<complex(double)> dtm;
    std::function<complex(double, double)> dtm2d;

    std::vector<double> ts, psx, psy, psz;
    std::vector<double> frequencies;
    std::vector<complex> dipole;
    std::vector<std::vector<complex>> dipole2d;
    std::vector<complex> hhg;



    void SetupTimeIntegrationVariables(double dt, double tmax);
    void SetupMomentumIntegrationVariables(double dp, double pmax, double pmin);
    void SetupMomentumIntegrationVariables(double dpx, double dpy, double pxmax, double pymax);
    void SetupField();
    void SetupFrequencyVariables(double df, double pmax, double pmin);

    double Action(double px, double py, double pz, int timstep) const;
    double Action(double px, double py, int timestep) const;
    double Action(double px, int timestep) const;
    void Execute1d();
    void Execute2d();

    void ComputeHHG1D();
    bool StoreField(const std::string& filename);
    bool StoreDipole(const std::string& filename);
    bool StoreHHG1D(const std::string& filename);
    

};
