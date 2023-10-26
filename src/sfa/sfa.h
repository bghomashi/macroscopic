#pragma once

#include <complex>
#include <vector>
#include <functional>
#include <memory>

using namespace std::complex_literals;
typedef std::complex<double> complex;

struct Pulse {
    typedef std::shared_ptr<Pulse> Ptr_t;
    virtual double A(double t) const = 0;
    virtual double E(double t) const = 0;
    virtual double alpha(double t) const = 0;
    virtual double beta(double t) const = 0;
    
    bool Store(const std::string& filename, const std::vector<double>& time) const;
    bool StoreA(const std::string& filename, const std::vector<double>& times) const;
    bool StoreAA(const std::string& filename, const std::vector<double>& times) const;
};

struct SFA {
    double Ip;

    Pulse::Ptr_t pulse;
    std::function<complex(double)> dtm;

    std::vector<double> trs, tis, ps;
    std::vector<double> frequencies;
    std::vector<complex> dipole;
    std::vector<complex> hhg;


    void SetupIonizationTimeIntegrationVariables(double dt, double tmax);
    void SetupRecombinationTimeIntegrationVariables(double dt, double tmax);
    void SetupMomentumIntegrationVariables(double dp, double pmax, double pmin);
    void SetupFrequencyVariables(double df, double pmax, double pmin);

    double Action(double p, double t) const;
    void Execute();

    void ComputeHHG();
    void ComputeMomentumTimeDistribution(const std::string& filename);
    void ComputeMomentumFreqDistribution(const std::string& filename);
};