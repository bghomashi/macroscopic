#include "sfa.h"
#include "integrate.h"
#include "maths.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

std::vector<double> blackman(int N) {
    std::vector<double> w(N);
    int M = (N % 2 == 0 ? N/2 : (N+1)/2);
    for (int n = 0; n <= M-1; n++) 
        w[n] = 0.42 - 0.5*cos(2*pi*n / (N-1)) + 0.08*cos(4*pi*n/(N-1));
    for (int n = M; n <= N-1; n++) 
        w[n] = 0.42 - 0.5*cos(2*pi*(N-1-n) / (N-1)) + 0.08*cos(4*pi*(N-1-n)/(N-1));
    return w;
}

void SFA::SetupTimeIntegrationVariables(double dt, double tmax) {
    int nt = tmax / dt + 1;
    ts.resize(nt);
    for (int i = 0; i < nt; i++)
        ts[i] = i*dt;
}
void SFA::SetupMomentumIntegrationVariables(double dp, double pmax, double pmin) {
    int np = (pmax - pmin) / dp + 1;
    psx.resize(np);
    for (int i = 0; i < np; i++)
        psx[i] = pmin + i*dp;
}
void SFA::SetupMomentumIntegrationVariables(double dpx, double dpy, double pxmax, double pymax){
    int npx = (pxmax) / dpx + 1;
    int npy = (pymax) / dpy + 1;
    int np = npx * npy;
    psx.resize(npx);
    for (int i = 0; i < npx; i++)
        psx[i] = i*dpx;
    psy.resize(npy);
    for (int i = 0; i < npy; i++)
        psy[i] = i*dpy;
}
void SFA::SetupFrequencyVariables(double df, double fmax, double fmin) {
    int nf = (fmax - fmin) / df + 1;
    frequencies.resize(nf);
    for (int i = 0; i < nf; i++)
        frequencies[i] = fmin + i*df;
}

void SFA::SetupField() {
    ///functions for computing integral of total field, and square of total field
    std::function<double(double)> Axtot = [this](double t) {
        double axtot = 0;
        for (auto& p : pulse)
            axtot += p->A(t).x;
        return axtot;
    };
    std::function<double(double)> Asqxtot = [this](double t) {
        double asqxtot = 0;
        for (auto& p : pulse)
            asqxtot += p->A(t).x * p->A(t).x;
        return asqxtot;
    };
    std::function<double(double)> Aytot = [this](double t) {
        double aytot = 0;
        for (auto& p : pulse)
            aytot += p->A(t).y;
        return aytot;
    };
    std::function<double(double)> Asqytot = [this](double t) {
        double asqytot = 0;
        for (auto& p : pulse)
            asqytot += p->A(t).y * p->A(t).y;
        return asqytot;
    };
    std::function<double(double)> Aztot = [this](double t) {
        double aztot = 0;
        for (auto& p : pulse)
            aztot += p->A(t).z;
        return aztot;
    };
    std::function<double(double)> Asqztot = [this](double t) {
        double asqztot = 0;
        for (auto& p : pulse)
            asqztot += p->A(t).z * p->A(t).z;
        return asqztot;
    };
    integralA.resize(ts.size());
    integralAsq.resize(ts.size());
    
    //produces time indexed vectors for integral of A and A^2
    std::vector<double> integralAx = TrapzInd(ts, Aztot);
    std::vector<double> integralAy = TrapzInd(ts, Aytot);
    std::vector<double> integralAz = TrapzInd(ts, Aztot);
    std::vector<double> integralAx2 = TrapzInd(ts, Asqxtot);
    std::vector<double> integralAy2 = TrapzInd(ts, Asqytot);
    std::vector<double> integralAz2 = TrapzInd(ts, Asqztot);
    for (int i = 0; i < ts.size(); i++) {
        integralA[i].x = integralAx[i];
        integralA[i].y = integralAy[i];
        integralA[i].z = integralAz[i];
        integralAsq[i].x = integralAx2[i];
        integralAsq[i].y = integralAy2[i];
        integralAsq[i].z = integralAz2[i];
    } 
    
    //produced time indexed vector for E and A
    Etot.resize(ts.size(), {0, 0, 0});
    Atot.resize(ts.size(), {0, 0, 0});
    for (int i = 0; i < ts.size(); i++) {
        for (int j = 0; j < pulse.size(); j++){
            Etot[i].x += pulse[j]->E(ts[i]).x;
            Etot[i].y += pulse[j]->E(ts[i]).y;
            Etot[i].z += pulse[j]->E(ts[i]).z;
            Atot[i].x += pulse[j]->A(ts[i]).x;
            Atot[i].y += pulse[j]->A(ts[i]).y;
            Atot[i].z += pulse[j]->A(ts[i]).z;
        }
    }
    
    
}

double SFA::Action(double px, double py, double pz, int timstep) const {
 ///action within time integral is integral[[p+A]^2, dt] + Ip*t
    double t = ts[timstep];
    double H0 = px*px + py*py + pz*pz + Ip;
    double pA = 2*px*integralA[timstep].x + py*integralA[timstep].y + pz*integralA[timstep].z;
    double AA = integralAsq[timstep].x + integralAsq[timstep].y + integralAsq[timstep].z;
    return .5*(H0 * t + pA + AA);

}
double SFA::Action(double px, double py, int timestep) const
{
    double t = ts[timestep];
    double H0 = px*px + py*py + Ip;
    double pA = 2*px*integralA[timestep].x + py*integralA[timestep].y;
    double AA = integralAsq[timestep].x + integralAsq[timestep].y;
    return .5*(H0 * t + pA + AA);
}
double SFA::Action(double px, int timestep) const {
    double t = ts[timestep];
    double H0 = px*px + Ip;
    double pA = 2*px*integralA[timestep].x;
    double AA = integralAsq[timestep].x;
    return .5*(H0 * t + pA + AA);
}
void SFA::Execute1d() {
    
    dipole.resize(ts.size(), 0);
    std::function<complex(double,int)> timeIntegrand = [this](double p, int t){
        return exp(1i*(Action(p, t))) *  Etot[t].x * dtm(p + Atot[t].x);
    };

    std::vector<complex> timeintegral;
    std::vector<complex> momentumIntegrand;
    momentumIntegrand.resize(psx.size(),0);
    timeintegral.resize(psx.size(),0);

    for (int i = 1; i < ts.size(); i++){
        for (int j = 0; j < psx.size(); j++){
            timeintegral[j] = timeintegral[j] + 0.5*(timeIntegrand(psx[j], i) + timeIntegrand(psx[j], i-1))*(ts[i] - ts[i-1]);
            momentumIntegrand[j] = timeintegral[j] * exp(-1i*Action(psx[j], i)) * std::conj(dtm(psx[j] + Atot[i].x));
        }
        //calculate the integral over momentum at the time step
        dipole[i] = Trapz(psx, momentumIntegrand);
    }
}

void SFA::Execute2d() {
    dipole2d.resize(ts.size(), std::vector<complex>(2, 0));
    std::function<complex(double,double,int)> timeIntegrandx = [this](double px, double py, int t){
        return exp(1i*Action(px, py, t)) *  Etot[t].x * dtm2d(px + Atot[t].x, py + Atot[t].y);
    };
    std::function<complex(double, double, int)> timeIntegrandy = [this](double px, double py, int t){
        return exp(1i*Action(px, py, t)) *  Etot[t].y * dtm2d(px + Atot[t].x, py + Atot[t].y);
    };

    std::function<std::vector<complex>(double, double, int)> timeIntegrand = [this](double px, double py, int t){
        std::vector<complex> timeIntegrand;
        timeIntegrand.resize(2);
        timeIntegrand[0] = exp(1i*Action(px, py, t)) *  Etot[t].x * dtm2d(px + Atot[t].x, py + Atot[t].y);
        timeIntegrand[1] = exp(1i*Action(px, py, t)) *  Etot[t].y * dtm2d(px + Atot[t].x, py + Atot[t].y);
        return timeIntegrand;
    };

    std::vector<std::vector<std::vector<complex>>> timeintegral;
    std::vector<std::vector<std::vector<complex>>> momentumIntegrand; //data structure looks like f[0](x)(y) is the x component and f[1](x)(y) is the y component
    momentumIntegrand.resize(2, std::vector<std::vector<complex>>(psx.size(), std::vector<complex>(psy.size(), 0)));
    timeintegral.resize(2, std::vector<std::vector<complex>>(psx.size(), std::vector<complex>(psy.size(), 0)));

    for (int i = 1; i < ts.size(); i++){ //for each time
        for(int j =0; j < 2; j++){ //for each component of the field
            for(int k = 0; k < psx.size(); k++){ //for each x momentum
                for(int l =0; l < psy.size(); l++){ //for each y momentum
                    //perform trapezoidal rule to get the next timestep
                    timeintegral[j][k][l] = timeintegral[j][k][l] + 0.5*(timeIntegrand(psx[k], psy[l], i)[j] + timeIntegrand(psx[k], psy[l], i-1)[j])*(ts[i] - ts[i-1]);
                    momentumIntegrand[j][k][l] = timeintegral[j][k][l] * exp(-1i*Action(psx[k], psy[l], i)) * std::conj(dtm2d(psx[j] + Atot[i].x, psy[k] + Atot[i].y));
                }
            }
        }
        for(int j = 0; j < 2; j++) //for each component of the field
            //calculate the integral over momentum at the time step
            dipole2d[i][j] = Trapz2D(psx, psy, momentumIntegrand[j]);
    }
    
}
void SFA::ComputeHHG1D() {
    //window dipole in time domain
    /*
    std::vector<double> window = blackman(dipole.size());
    for (int i = 0; i < dipole.size(); i++) {
        dipole[i] *= window[i];
    }
    */
    //compute hhg
    hhg.resize(frequencies.size(), 0);
    for (int i = 0; i < frequencies.size(); i++) {
        for (int j = 0; j < ts.size(); j++) {
            hhg[i] += dipole[j] * exp(-1i * frequencies[i] * ts[j]);
        }
    } 
}

bool SFA::StoreHHG1D(const std::string& filename){
    int nf = frequencies.size();
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int iff = 0; iff < nf; iff++) {
        file << frequencies[iff] << "\t"
             << std::real(hhg[iff]) << "\t"
             << std::imag(hhg[iff]) << "\n";

    }
    return true;
}









bool Pulse::Store(const std::string& filename, const std::vector<double>& times) const {
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int i = 0; i < times.size(); i++) {
        file << times[i] << "\t"
             << A(times[i]).x << "\t"
             << A(times[i]).y << "\t"
             << A(times[i]).z << "\t"
             << E(times[i]).x << "\t"
             << E(times[i]).y << "\t"    
             << E(times[i]).z << "\n";
    }

    return true;
}


