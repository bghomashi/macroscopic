#include "sfa.h"
#include "sin2pulse.h"
#include "maths.h"
#include "../maths/constants.h"


// pulse params
double E0 = sqrt(1e14/3.51e16);
double w0 = 0.057;
double N = 10;
double CEP = 0;

// sfa params
double dt = 0.1;
double tmax = (2*pi/w0) * N;
double dp = 0.01;
double pmin = -3, pmax = 3;
double dff = 0.01;
double ffmin = 0.*w0, ffmax = 40.*w0;


int main() {
    
    SFA sfa;
    
    sfa.Ip = 0.5;
    sfa.pulse.push_back(std::make_shared<Sin2Pulse>(CEP, w0, N, E0));
    
    sfa.dtm = [&sfa](double p) {
        double K = sqrt(2.*sfa.Ip);
        return -2.i*sqrt(K)*p / (p*p + K*K) / (p*p + K*K);
    };
    
    sfa.SetupTimeIntegrationVariables(dt, tmax);
    sfa.SetupMomentumIntegrationVariables(dp, pmax, pmin);
    sfa.SetupFrequencyVariables(dff, ffmax, ffmin);

    // sfa.StoreAction("action.dat");

    // sfa.pulse->Store("pulse.dat", sfa.trs);
    // sfa.pulse->StoreA("pulseA.dat", sfa.trs);
    // sfa.pulse->StoreAA("pulseAA.dat", sfa.trs);
    
    // sfa.ComputeMomentumTimeDistribution("momentum.dat");
    // sfa.ComputeMomentumFreqDistribution("momentum_freq.dat");
    sfa.SetupField();
    
    sfa.Execute1d();
    
    //sfa.StoreDipole("dipole.out");
    // sfa.StoreDipoleWindowed("dipole_windowed.out");

    sfa.ComputeHHG1D();
    sfa.StoreHHG1D("hhg_At.out");
    sfa.StoreDipole("dipole.out");

    
    return 0;
}
