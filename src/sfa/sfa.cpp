#include "sfa.h"
#include "integrate.h"
#include "maths.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <chrono>

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
void SFA::SetupMomentumIntegrationVariables(double dpx, double pxmax, double pxmin, double dpy, double pymax, double pymin) {
    int npx = (pxmax - pxmin) / dpx + 1;
    int npy = (pymax - pymin) / dpy + 1;
    psx.resize(npx);
    psy.resize(npy);
    for (int i = 0; i < npx; i++)
        psx[i] = pxmin + i*dpx;
    for (int i = 0; i < npy; i++)
        psy[i] = pymin + i*dpy;
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
    std::function<double(double)> Asqxtot = [this, Axtot](double t) {
        return (Axtot(t) * Axtot(t));
    };
    std::function<double(double)> Aytot = [this](double t) {
        double aytot = 0;
        for (auto& p : pulse)
            aytot += p->A(t).y;
        return aytot;
    };
    std::function<double(double)> Asqytot = [this, Aytot](double t) {
        return (Aytot(t) * Aytot(t));
    };
    std::function<double(double)> Aztot = [this](double t) {
        double aztot = 0;
        for (auto& p : pulse)
            aztot += p->A(t).z;
        return aztot;
    };
    std::function<double(double)> Asqztot = [this, Aztot](double t) {
        return (Aztot(t) * Aztot(t));
    };
    integralA.resize(ts.size());
    integralAsq.resize(ts.size());
    
    //produces time indexed vectors for integral of A and A^2
    std::vector<double> integralAx = TrapzInd(ts, Axtot);
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
double SFA::Action(double px, int t, int tp) const {
    return Action(px, t) - Action(px, tp);
}
double SFA::Action(double px, double py, int t, int tp) const {
    return Action(px, py, t) - Action(px, py, tp);
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


void SFA::SaddlePoint1d()
{
    dipole.resize(ts.size(), 0);
    double spm; //saddle point momentum
    complex prev = 0;// previous time integrand
    complex next = 0;// next time integrand

    std::function<double(double t, double tp)> saddleMomentum = [this](double t, double tp){
        return (integralA[t].x - integralA[tp].x)/(ts[t] - ts[tp]);
    }; // saddle point momentum

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

    for (int i = 1; i < ts.size(); i++){ //t
        for (int j = 1; j < i; j++){ //t'
            spm = saddleMomentum(i, j); 
            next = std::conj(dtm(spm + Atot[i].x)) * exp(-1i*Action(spm, i, j)) * Etot[j].x * dtm(spm + Atot[j].x);
            dipole[i] += .5 * (prev  + next) * (ts[j] - ts[j-1]); // trapezoidal rule
            prev = next;
        }
        //-------------------for tracking progress of the calculation-------------------
        std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> timeSpan = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
        double iterationsPerSecond = i / timeSpan.count();
        double estimatedTime = (ts.size() - i) / iterationsPerSecond;
        int minutes = estimatedTime / 60;
        int hours = minutes / 60;
        int estimatedTimeSeconds = (int)estimatedTime % 60;
        int estimatedTimeMinutes = minutes % 60;
        int estimatedTimeHours = hours;

        std::cout << std::setprecision(2) << "its[" << i << " / " << ts.size() << "] its/s[" << iterationsPerSecond
                  << "] remaining[" << estimatedTimeHours << "h " << estimatedTimeMinutes << "m " << estimatedTimeSeconds << "s]"  << "\r" << std::flush;
    }
}

void SFA::SaddlePoint2d() {
    dipole2d.resize(ts.size(), std::vector<complex>(2, 0));
    double spmx, spmy; //saddle point momentum
    complex prevx = 0;// previous time integrand
    complex prevy = 0;// previous time integrand
    complex nextx = 0;// next time integrand
    complex nexty = 0;// next time integrand
    complex EdotD; 

    std::function<double(double t, double tp)> saddleMomentumx = [this](double t, double tp) {
        return (integralA[t].x - integralA[tp].x) / (ts[t] - ts[tp]);
    }; // saddle point momentum
    std::function<double(double t, double tp)> saddleMomentumy = [this](double t, double tp) {
        return (integralA[t].y - integralA[tp].y) / (ts[t] - ts[tp]);
    }; // saddle point momentum

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    
    for (int i = 1; i < ts.size(); i++) { //t
        for (int j = 1; j < i; j++) { //t'
            spmx = saddleMomentumx(i, j);
            spmy = saddleMomentumy(i, j);
            EdotD = Etot[j].x * dtm2d(spmx + Atot[j].x, spmy + Atot[j].y)[0] + Etot[j].y * dtm2d(spmx + Atot[j].x, spmy + Atot[j].y)[1];
            nextx = std::conj(dtm2d(spmx + Atot[i].x, spmy + Atot[i].y)[0]) * exp(-1i*Action(spmx, spmy, i, j)) * EdotD;
            nexty = std::conj(dtm2d(spmx + Atot[i].x, spmy + Atot[i].y)[1]) * exp(-1i*Action(spmx, spmy, i, j)) * EdotD;
            dipole2d[i][0] += .5 * (prevx + nextx) * (ts[j] - ts[j - 1]); // trapezoidal rule
            dipole2d[i][1] += .5 * (prevy + nexty) * (ts[j] - ts[j - 1]); // trapezoidal rule
            prevx = nextx;
            prevy = nexty;
        }
        
        //-------------------for tracking progress of the calculation-------------------
        std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> timeSpan = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
        double iterationsPerSecond = i / timeSpan.count();
        double estimatedTime = (ts.size() - i) / iterationsPerSecond;
        int minutes = estimatedTime / 60;
        int hours = minutes / 60;
        int estimatedTimeSeconds = (int)estimatedTime % 60;
        int estimatedTimeMinutes = minutes % 60;
        int estimatedTimeHours = hours;

        std::cout << std::setprecision(2) << "its[" << i << " / " << ts.size() << "] its/s[" << iterationsPerSecond
                << "] remaining[" << estimatedTimeHours << "h " << estimatedTimeMinutes << "m " << estimatedTimeSeconds << "s]"  << "\r" << std::flush; 
    }
}

void SFA::Execute2d() {

    
    dipole2d.resize(ts.size(), std::vector<complex>(2, 0));
    std::function<complex(double, double, int)> timeIntegrand = [this](double px, double py, int t){
        return exp(1i*Action(px, py, t)) * (Etot[t].x * dtm2d(px + Atot[t].x, py + Atot[t].y)[0] + Etot[t].y * dtm2d(px + Atot[t].x, py + Atot[t].y)[1]);
    };
    std::vector<std::vector<complex>> timeintegral;
    std::vector<std::vector<std::vector<complex>>> momentumIntegrand; //data structure looks like f[0](x)(y) is the x component and f[1](x)(y) is the y component
    momentumIntegrand.resize(2, std::vector<std::vector<complex>>(psx.size(), std::vector<complex>(psy.size(), 0)));
    //resize time integral to be an array that is psx.size by psy.size
    timeintegral.resize(psx.size(), std::vector<complex>(psy.size(), 0));
    

    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

    for (int i = 1; i < ts.size(); i++){ //for each time
        for(int j =0; j < 2; j++){ //for each component of the field
            for(int k = 0; k < psx.size(); k++){ //for each x momentum
                for(int l =0; l < psy.size(); l++){ //for each y momentum
                    //perform trapezoidal rule to get the next timestep
                    timeintegral[k][l] = timeintegral[k][l] + 0.5*(timeIntegrand(psx[k], psy[l], i) + timeIntegrand(psx[k], psy[l], i-1))*(ts[i] - ts[i-1]);
                    momentumIntegrand[j][k][l] = timeintegral[k][l] * exp(-1i*Action(psx[k], psy[l], i)) * std::conj(dtm2d(psx[j] + Atot[i].x, psy[k] + Atot[i].y)[j]);
                }
            }
            //calculate the integral over momentum at the time step
            dipole2d[i][j] = Trapz2D(psx, psy, momentumIntegrand[j]);
        }
        //-------------------for tracking progress of the calculation-------------------
        std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> timeSpan = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
        double iterationsPerSecond = i / timeSpan.count();
        double estimatedTime = (ts.size() - i) / iterationsPerSecond;
        int minutes = estimatedTime / 60;
        int hours = minutes / 60;
        int estimatedTimeSeconds = (int)estimatedTime % 60;
        int estimatedTimeMinutes = minutes % 60;
        int estimatedTimeHours = hours;

        std::cout << std::setprecision(2) << "its[" << i << " / " << ts.size() << "] its/s[" << iterationsPerSecond
                  << "] remaining[" << estimatedTimeHours << "h " << estimatedTimeMinutes << "m " << estimatedTimeSeconds << "s]"  << "\r" << std::flush;

    }
    std::cout << std::endl;
    
}

void SFA::ComputeHHG1D() {
    //window dipole in time domain

    std::vector<double> window = blackman(dipole.size());
    for (int i = 0; i < dipole.size(); i++) {
        dipole[i] *= window[i];
    }
    //compute hhg
    hhg.resize(frequencies.size(), 0);
    for (int i = 0; i < frequencies.size(); i++) {
        for (int j = 0; j < ts.size(); j++) {
            hhg[i] += dipole[j] * exp(-1i * frequencies[i] * ts[j]);
        }
    } 
}

void SFA::ComputeHHG2D() {
    //window dipole in time domain

    std::vector<double> window = blackman(dipole2d.size());
    for (int i = 0; i < dipole2d.size(); i++) {
        dipole2d[i][0] *= window[i];
        dipole2d[i][1] *= window[i];
    }
    hhg2d.resize(frequencies.size(), std::vector<complex>(2, 0));
    for (int i = 0; i < frequencies.size(); i++) {
        for (int j = 0; j < ts.size(); j++) {
            hhg2d[i][0] += dipole2d[j][0] * exp(-1i * frequencies[i] * ts[j]);
            hhg2d[i][1] += dipole2d[j][1] * exp(-1i * frequencies[i] * ts[j]);
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

bool SFA::StoreHHG2D(const std::string& filename)
{

    int nf = frequencies.size();
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;
    for (int iff = 0; iff < nf; iff++) {
        file << frequencies[iff] << "\t"
             << std::real(hhg2d[iff][0]) << "\t"
             << std::imag(hhg2d[iff][0]) << "\t"
             << std::real(hhg2d[iff][1]) << "\t"
             << std::imag(hhg2d[iff][1]) << "\n";
    }

}

bool SFA::StoreDipole(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;
    for(int i = 0; i < dipole.size(); i++){
        file << ts[i] << "\t"
             << std::real(dipole[i]) << "\t"
             << std::imag(dipole[i]) << "\n";
    }
}

bool SFA::StoreDipole2d(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;
    for (int i = 0; i < dipole2d.size(); i++) {
        file << ts[i] << "\t"
             << std::real(dipole2d[i][0]) << "\t"
             << std::imag(dipole2d[i][0]) << "\t"
             << std::real(dipole2d[i][1]) << "\t"
             << std::imag(dipole2d[i][1]) << "\n";
    }
}

bool SFA::StoreDTM2d(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;
    for (int i = 0; i < psx.size(); i++) {
        for (int j = 0; j < psy.size(); j++) {
            file << psx[i] << "\t"
                 << psy[j] << "\t"
                 << std::real(dtm2d(psx[i], psy[j])[0]) << "\t"
                 << std::imag(dtm2d(psx[i], psy[j])[0]) << "\t"
                 << std::real(dtm2d(psx[i], psy[j])[1]) << "\t"
                 << std::imag(dtm2d(psx[i], psy[j])[1]) << "\n";
        }
    }
}

bool SFA::StoreField(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;
    for (int i = 0; i < ts.size(); i++) {
        file << ts[i] << "\t"
             << Etot[i].x << "\t"
             << Etot[i].y << "\t"
             << Etot[i].z << "\t"
             << Atot[i].x << "\t"
             << Atot[i].y << "\t"
             << Atot[i].z << "\t"
             << integralA[i].x << "\t"
             << integralA[i].y << "\t"
             << integralA[i].z << "\t"
             << integralAsq[i].x << "\t"
             << integralAsq[i].y << "\t"
             << integralAsq[i].z << "\n";


    }
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





