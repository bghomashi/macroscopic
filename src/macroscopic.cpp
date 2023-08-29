#include "macroscopic.h"
#include "utility/profiler.h"
#include <iostream>

cvector Macroscopic::Spectrum(double td, double pd) const {
    cvector spectrum(_microscopic_data.Frequencies().size(), 0);    //initialize spectrum to 0
    cvector microspec(_microscopic_data.Frequencies().size(), 0);   // storage for microscopic...
    auto& frequencies = _microscopic_data.Frequencies();

    int count = 0;
    for (auto& cell : _points._cells) {
        // skip intensities too weak
        if (cell.intensity / _laser.PeakI0() < 0.032) {
            count++;
            continue;
        }
        auto p = cell.pos;
        cvector microspec = _microscopic_data.Lerp(cell.intensity);
        for (int j = 0; j < frequencies.size(); j++) {
            complex radiation = E(p, td, pd, frequencies[j], microspec[j]);
            spectrum[j] += cell.density * radiation;
        }
    }
    std::cout << "Count below cutoff: " << count << std::endl;
    return spectrum;
}
complex Macroscopic::E( const point3& rj,                           // dipole position
                        const double theta_d, const double phi_d,   // detector location
                        double w, complex aj) const {               // frequency and radiation coefficient
    double arg =    rj.x * cos(theta_d)*sin(phi_d) + 
                    rj.y * sin(theta_d)*sin(phi_d) + 
                    rj.z * (cos(theta_d) - 1);

    double n = 2*std::round((w / _laser.Freq() - 1) / 2) + 1;
    double phase = _laser.Phase(rj);


    return aj * exp(-1.i * (w / C) * arg) * exp(-1.i*n*phase);
}



void Macroscopic::Initialize(Laser& laser, points samples, const std::string& filename) {
    
    ProfilerPush();

    // load microscopic data
    _microscopic_data = SpectrumMatrix::Load(filename);
    // setup laser
    _laser = laser; 
    // setup gas jet
    _points = samples;
    _points.CalcCells(_laser);

    ProfilerPop();
}

dvector Macroscopic::Frequencies() const {
    return _microscopic_data.Frequencies();
}
cvector Macroscopic::MicroSpectrum(const point3& position) const {
    double I = _laser.IntensityAt(position);

    return _microscopic_data.Lerp(I);
}