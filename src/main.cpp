#include <iostream>
#include <fstream>
#include <iomanip>
#include "macroscopic.h"
#include "utility/profiler.h"

double wavelength_nm = 800;
double beam_waist_um = 10;
double peak_I0_wcm2 = 1e14;
double gas_radius = 500;
double gas_length = 2;
size_t number_of_cells = 1e5;
double jetsig_um = 50;
double density_cm3 = 1e18;
std::string filename = "hhg_vs_int.txt";
std::vector<std::pair<double, double>> detectors = {{0, 0}};
std::vector<cvector> spectrums(detectors.size());

int main() {
    Macroscopic macroscopic;
    macroscopic.Initialize(
        wavelength_nm, beam_waist_um, peak_I0_wcm2, 
        gas_radius, gas_length, number_of_cells,
        jetsig_um, density_cm3,
        filename
    );
    auto H = macroscopic.Frequencies();
    for (int i = 0; i < detectors.size(); i++) {
        auto& d = detectors[i];
        double td = d.first;
        double pd = d.second;

        Profile::Push("Spectrum");
        spectrums[i] = macroscopic.Spectrum(td, pd);
        Profile::Pop("Spectrum");
    }

    std::ofstream file("data.txt");
    file << std::setprecision(8) << std::scientific;
    // for (int i = 0; i < detectors.size(); i++)
    //     file << i << "\t";
    // file << "\n";

    for (int i = 0; i < H.size(); i++) {
        file << H[i] << "\t";
        for (int j = 0; j < detectors.size(); j++) {
            file << abs(spectrums[j][i])*abs(spectrums[j][i]) << "\t";
        }
        file << "\n";
    }


    Profile::Print();

    return 0;
}