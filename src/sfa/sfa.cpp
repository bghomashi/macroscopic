#include "sfa.h"
#include "integrate.h"
#include "maths.h"
#include <iostream>
#include <fstream>
#include <iomanip>

std::vector<double> blackman(int N) {
    std::vector<double> w(N);
    int M = (N % 2 == 0 ? N/2 : (N+1)/2);
    for (int n = 0; n <= M-1; n++) 
        w[n] = 0.42 - 0.5*cos(2*pi*n / (N-1)) + 0.08*cos(4*pi*n/(N-1));
    for (int n = M; n <= N-1; n++) 
        w[n] = 0.42 - 0.5*cos(2*pi*(N-1-n) / (N-1)) + 0.08*cos(4*pi*(N-1-n)/(N-1));
    return w;
}

void SFA::SetupIonizationTimeIntegrationVariables(double dt, double tmax) {
    int nt = tmax / dt + 1;
    tis.resize(nt);
    for (int i = 0; i < nt; i++)
        tis[i] = i*dt;
}
void SFA::SetupRecombinationTimeIntegrationVariables(double dt, double tmax) {
    int nt = tmax / dt + 1;
    trs.resize(nt);
    for (int i = 0; i < nt; i++)
        trs[i] = i*dt;
}
void SFA::SetupMomentumIntegrationVariables(double dp, double pmax, double pmin) {
    int np = (pmax - pmin) / dp + 1;
    ps.resize(np);
    for (int i = 0; i < np; i++)
        ps[i] = pmin + i*dp;
}
void SFA::SetupFrequencyVariables(double df, double fmax, double fmin) {
    int nf = (fmax - fmin) / df + 1;
    frequencies.resize(nf);
    for (int i = 0; i < nf; i++)
        frequencies[i] = fmin + i*df;
}
double SFA::Action(double p, double t) const {
    double H0 = 0.5 * p * p + Ip;
    double pA = p * pulse->alpha(t);
    double AA = 0.5 * pulse->beta(t);

    return H0 * t + pA + AA;
}
void SFA::Execute() {
    int64_t ntr = trs.size();
    int64_t nti = tis.size();
    int64_t np = ps.size();
    
    this->dipole.resize(ntr, 0);


    // function to compute *ionization time* dependent integrand
    auto Gamma = [this] (double p, double ti) {
        return  exp(1i*Action(p, ti)) *
                pulse->E(ti) *
                dtm(p + pulse->A(ti));
    };

    // returns the integral up to each time step for a given momentum p
    auto compute = [this, Gamma] (double p) {
        int ntr = trs.size();
        int nti = tis.size();
        std::vector<complex> GG(ntr);
        int iti = 1;
        GG[0] = 0;
        for (int itr = 1; itr < ntr; itr++) {                           // for each recombination time, tr
            double tr = trs[itr];                                       // this time
            GG[itr] = GG[itr-1];                                        // start with previous value
            for (; iti < nti && tis[iti] <= tr; iti++)                  // for each ionization time BEFORE tr
                GG[itr] += 0.5 * (Gamma(p, tis[iti]) + Gamma(p, tis[iti-1]))*(tis[iti] - tis[iti-1]);   // trap
        }
        return GG;
    };
    
    std::vector<complex> GG0(ntr), GG1(ntr);
    GG0 = compute(ps[0]);
    for (int ip = 1; ip < np; ip++) {
        std::cout << "step: " << ip << "/" << np << std::flush;

        double p0 = ps[ip - 1];
        double p1 = ps[ip];

        GG1 = compute(ps[ip]);

        for (int it = 0; it < ntr; it++) {
            double tr = trs[it];
            complex W0 = 0;
            complex W1 = 0;
            W0 = std::conj(dtm(p0 + pulse->A(tr))) * exp(-1i*Action(p0, tr)) *  GG0[it];
            W1 = std::conj(dtm(p1 + pulse->A(tr))) * exp(-1i*Action(p1, tr)) *  GG1[it];
            
            // if (p0 + pulse->A(tr) != 0)
            //     W0 = std::conj(dtm(p0 + pulse->A(tr))) * exp(-1i*Action(p0, tr)) *  GG0[it] * (pulse->A(tr) / (p0 + pulse->A(tr)));
            // if (p1 + pulse->A(tr) != 0)
            //     W1 = std::conj(dtm(p1 + pulse->A(tr))) * exp(-1i*Action(p1, tr)) *  GG1[it] * (pulse->A(tr) / (p1 + pulse->A(tr)));
            //complex W0 = std::conj(dtm(p0 + pulse->A(tr))); //* GG0[it];
            //complex W1 = std::conj(dtm(p1 + pulse->A(tr))); //* GG1[it];
            //complex W0 = std::conj(dtm(p0)) /* exp(-1i*Action(p0, tr))*/ *  GG0[it];
            //complex W1 = std::conj(dtm(p1)) /* exp(-1i*Action(p1, tr))*/ *  GG1[it];
            //complex W0 = GG0[it];
            //complex W1 = GG1[it];

            dipole[it] += 0.5*(W1 + W0)*(p1 - p0);

            // displays progress count
            if ((it == ntr / 3) || (it == 2 * ntr / 3) || (it == ntr -1))
                std::cout << "." << std::flush;
        }
        
        std::swap(GG0, GG1);
        
        std::cout << "Complete." << std::endl;
    }
    for (int it = 0; it < ntr; it++) {
        dipole[it] *= -1.i;
    }
    for (int it = 0; it < ntr; it++) {
        dipole[it] += std::conj(dipole[it]);
    }
}
bool SFA::StoreDipole(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open())
        return false;

    int nt = trs.size();
    file << std::setprecision(8) << std::scientific;

    for (int it = 0; it < nt; it++)
        file << trs[it] << "\t" 
             << std::real(dipole[it]) << "\t" 
             << std::imag(dipole[it]) << std::endl;

    return true;
}
bool SFA::StoreDipoleWindowed(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open())
        return false;

    int nt = trs.size();
    auto window = blackman(nt);
    file << std::setprecision(8) << std::scientific;

    for (int it = 0; it < nt; it++)
        file << trs[it] << "\t" 
             << std::real(window[it]*dipole[it]) << "\t" 
             << std::imag(window[it]*dipole[it]) << std::endl;

    return true;
}

void SFA::ComputeHHG() {
    int nt = trs.size();
    int nf = frequencies.size();
    hhg.resize(nf);

    auto window = blackman(nt);

    for (int iff = 0; iff < nf; iff++) {
        double f = frequencies[iff];
        hhg[iff] = 0;
        for (int it = 0; it < nt; it++) {
            double t = trs[it];

            hhg[iff] += exp(-1i*f*t) * dipole[it];
        }
        hhg[iff] *= -f*f;
    }
}
bool SFA::StoreHHG(const std::string& filename) const {
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
             << A(times[i]) << "\t"
             << E(times[i]) << "\t"
             << alpha(times[i]) << "\t"
             << beta(times[i]) << "\n";
    }

    return true;
}
bool Pulse::StoreA(const std::string& filename, const std::vector<double>& times) const {
    std::vector<double> AA(times.size());
    for (int it = 0; it < times.size(); it++) {
        AA[it] = Trapz<double>(std::vector<double>(times.begin(), times.begin() + it), [this](double t){
            return A(t);
        });
    }

    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int i = 0; i < times.size(); i++) {
        file << times[i] << "\t"
             << AA[i] << "\n";
    }

    return true;
}
bool Pulse::StoreAA(const std::string& filename, const std::vector<double>& times) const {
    std::vector<double> AA(times.size());
    for (int it = 0; it < times.size(); it++) {
        AA[it] = Trapz<double>(std::vector<double>(times.begin(), times.begin() + it), [this](double t){
            return A(t)*A(t);
        });
    }

    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int i = 0; i < times.size(); i++) {
        file << times[i] << "\t"
             << AA[i] << "\n";
    }

    return true;
}
bool SFA::StoreAction(const std::string& filename) const {
    int64_t ntr = trs.size();
    int64_t np = ps.size();
    
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;


    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ntr; j++) {
            double S = Action(ps[i], trs[j]);
            file << ps[i] << "\t"
                 << trs[j] << "\t"
                 << S << "\n";
        }
        file << "\n";
    }

    return false;
}

void SFA::ComputeMomentumTimeDistribution(const std::string& filename) {
    int64_t ntr = trs.size();
    int64_t nti = tis.size();
    int64_t np = ps.size();
    
    std::vector<complex> mom(ntr*np, 0);

    // function to compute *ionization time* dependent integrand
    auto Gamma = [this] (double p, double ti) {
        return  exp(1i*Action(p, ti))*
                pulse->E(ti)*
                dtm(p + pulse->A(ti));
    };

    // returns the integral up to each time step for a given momentum p
    auto compute = [this, Gamma] (double p) {
        int ntr = trs.size();
        int nti = tis.size();
        std::vector<complex> GG(ntr);
        int iti = 1;
        GG[0] = 0;
        for (int itr = 1; itr < ntr; itr++) {                           // for each recombination time, tr
            double tr = trs[itr];                                       // this time
            GG[itr] = GG[itr-1];                                        // start with previous value
            for (; iti < nti && tis[iti] <= tr; iti++)                  // for each ionization time BEFORE tr
                GG[itr] += 0.5*(Gamma(p, tis[iti]) + Gamma(p, tis[iti-1]))*(tis[iti] - tis[iti-1]);   // trap
        }
        return GG;
    };
    
    std::vector<complex> GG1(ntr);
    for (int ip = 0; ip < np; ip++) {
        std::cout << "step: " << ip << "/" << np << std::flush;

        double p = ps[ip];

        GG1 = compute(ps[ip]);

        for (int it = 0; it < ntr; it++) {
            double tr = trs[it];
            complex W = std::conj(dtm(p + pulse->A(tr))) * exp(-1i*Action(p, tr)) *  GG1[it] / (p + pulse->A(tr));
            mom[ip*ntr + it] = W;
        }
                
        std::cout << "Complete." << std::endl;
    }




    std::ofstream file(filename);
    file << std::setprecision(8) << std::scientific;

    for (int i = 0; i < np; i++) {
        for (int j = 0; j < ntr; j++) {
            file << ps[i] /*- pulse->A(trs[j])*/ << "\t"
                 << trs[j] << "\t"
                 << std::real(mom[i*ntr + j]) << "\t"
                 << std::imag(mom[i*ntr + j]) << "\n";
        }
        file << "\n";
    }
}
void SFA::ComputeMomentumFreqDistribution(const std::string& filename) {
    int64_t ntr = trs.size();
    int64_t nti = tis.size();
    int64_t np = ps.size();
    int nf = frequencies.size();
    
    std::vector<complex> mom(ntr*np, 0);
    std::vector<complex> mom_w(nf*np, 0);

    // function to compute *ionization time* dependent integrand
    auto Gamma = [this] (double p, double ti) {
        return  exp(1i*Action(p, ti))*
                pulse->E(ti)*
                dtm(p + pulse->A(ti));
    };

    // returns the integral up to each time step for a given momentum p
    auto compute = [this, Gamma] (double p) {
        int ntr = trs.size();
        int nti = tis.size();
        std::vector<complex> GG(ntr);
        int iti = 1;
        GG[0] = 0;
        for (int itr = 1; itr < ntr; itr++) {                           // for each recombination time, tr
            double tr = trs[itr];                                       // this time
            GG[itr] = GG[itr-1];                                        // start with previous value
            for (; iti < nti && tis[iti] <= tr; iti++)                  // for each ionization time BEFORE tr
                GG[itr] += 0.5*(Gamma(p, tis[iti]) + Gamma(p, tis[iti-1]))*(tis[iti] - tis[iti-1]);   // trap
        }
        return GG;
    };
    
    std::vector<complex> GG1(ntr);
    for (int ip = 0; ip < np; ip++) {
        std::cout << "step: " << ip << "/" << np << std::flush;

        double p = ps[ip];

        GG1 = compute(ps[ip]);

        for (int it = 0; it < ntr; it++) {
            double tr = trs[it];
            complex W = std::conj(dtm(p + pulse->A(tr))) * exp(-1i*Action(p, tr)) * GG1[it] / (p + pulse->A(tr));
            mom[ip*ntr + it] = W;
        }
                
        std::cout << "Complete." << std::endl;
    }


    auto window = blackman(ntr);

    for (int iff = 0; iff < nf; iff++) {
        double f = frequencies[iff];
        for (int ip = 0; ip < np; ip++) {
            double p = ps[ip];
            mom_w[ip*nf + iff] = 0;
            for (int it = 0; it < ntr; it++) {
                double t = trs[it];
                
                if (p - pulse->A(t) < ps[0] || p - pulse->A(t) > ps.back())
                    continue;

                int index = std::upper_bound(ps.begin(), ps.end(), p - pulse->A(t)) - ps.begin();

                mom_w[ip*nf + iff] += exp(-1i*f*t) * 2. * window[it] * mom[ip*ntr + it] * (trs[1] - trs[0]);
            }
        }
    }



    std::ofstream file(filename);
    file << std::setprecision(8) << std::scientific;

    for (int j = 0; j < nf; j++) {
        for (int i = 0; i < np; i++) {
            file << frequencies[j] << "\t"
                 << ps[i] << "\t"
                 << std::real(mom_w[i*nf + j]) << "\t"
                 << std::imag(mom_w[i*nf + j]) << "\n";
        }
        file << "\n";
    }
}