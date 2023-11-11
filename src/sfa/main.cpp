#include "sfa.h"
#include "sin2pulse.h"
#include "maths.h"
#include "../maths/constants.h"
#include <iostream>

// pulse params
double E0 = sqrt(1e14/3.51e16);
double w0 = 0.057;
double N = 10;
double CEP = 0;

// sfa params
double dt = 0.05;
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

    sfa.dtm2d = [&sfa](double px, double py) {
        double K = sqrt(2.*sfa.Ip);
        return -2.i*sqrt(K *(px*px + py*py)) / (px*px + py*py + K*K) / (px*px + py*py + K*K);
    };
    
    sfa.SetupTimeIntegrationVariables(dt, tmax);
    sfa.SetupMomentumIntegrationVariables(dp, pmax, pmin);
    sfa.SetupFrequencyVariables(dff, ffmax, ffmin);
    sfa.SetupMomentumIntegrationVariables(dp, pmax, pmin, dp, pmax, pmin);
    sfa.SetupField();
    
    sfa.SaddlePoint2d();
    //sfa.SaddlePoint1d();
    //sfa.Execute1d();
    //sfa.Execute2d();

    //sfa.ComputeHHG1D();
    //sfa.StoreHHG1D("hhg_At.out");
    
    //sfa.StoreDipole2d("dipole.out");
    sfa.ComputeHHG2D();
    sfa.StoreHHG2D("hhg.out");
    std::cout << "Done." << std::endl;
    
    return 0;
}
