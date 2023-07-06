#include "macroscopic.h"
#include "utility/profiler.h"
#include <iostream>

cvector Macroscopic::Spectrum(double td, double pd) const {
    cvector spectrum(_microscopic_data.Frequencies().size(), 0);    //initialize spectrum to 0
    cvector microspec(_microscopic_data.Frequencies().size(), 0);   // storage for microscopic...
    auto& frequencies = _microscopic_data.Frequencies();

    int count = 0;
    for (auto& cell : _gas_jet._cells) {
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
    double rho_j = sqrt(rj.x*rj.x + rj.y*rj.y);         // location of dipole in laser
    double phase = _laser.Phase(rho_j, rj.z);


    return aj * exp(-1.i * (w / C) * arg) * exp(1.i*n*phase);
}



void Macroscopic::Initialize(
    double wavelength_nm, double beam_waist_um, double peakI0_wcm2, 
    double gas_radius_um, double gas_length_um, size_t number_of_cells,
    double jetsig_um, double density_cm3,
    const std::string& filename) {
    
    ProfilerPush();

    // convert to atomic units
    _jetsig_um = jetsig_um;
    _jetsig = umToAU * jetsig_um;
    _density_cm3 = density_cm3;
    _density = density_cm3 / cmToAU / cmToAU / cmToAU;
    
    // load microscopic data
    _microscopic_data = SpectrumMatrix::Load(filename);
    // setup laser
    _laser = Laser(peakI0_wcm2, beam_waist_um, wavelength_nm);
    // setup gas jet
    _gas_jet = CylindricalGasJet(_density, _jetsig, gas_radius_um*umToAU, gas_length_um*umToAU, number_of_cells);
    _gas_jet.SampleCylinder(_laser);
    // each point
    for (auto& cell : _gas_jet._cells) {
        double las_r = sqrt(cell.pos.x*cell.pos.x + cell.pos.y*cell.pos.y);
        double las_z = cell.pos.z;
        std::cout << cell.phase << std::endl;
    }

    ProfilerPop();
}
dvector Macroscopic::Frequencies() const {
    return _microscopic_data.Frequencies();
}
cvector Macroscopic::MicroSpectrum(const point3& position) const {
    double r = sqrt(position.x*position.x + position.y*position.y);
    double z = position.z;
    double I = _laser.IntensityAt(r, z);

    return _microscopic_data.Lerp(I);
}