#include <iostream>
#include <fstream>
#include <iomanip>
#include "macroscopic.h"
#include "utility/profiler.h"

#include "utility/json.hpp"

double peak_I0_wcm2 = 3.8e14;
double wavelength_nm = 800;
double beam_waist_um = 30;
double gas_radius = 500;
double gas_length = 6*beam_waist_um;
size_t number_of_cells = 1e5;
double jetsig_um = 800;
double density_cm3 = 1;
std::string filename = "he_hhg_vs_int.in";
std::vector<std::pair<double, double>> detectors = {{0,0}};
std::vector<cvector> spectrums(detectors.size());

int main() {
    // HOW TO READ JSON
    nlohmann::json input;                           // create json-object
    std::ifstream in_file_stream("input.json");     // open the (text) file
    in_file_stream >> input;                        // read the file into (>>) the json-object

    // HOW TO USE JSON
    std::cout << input["number"] << std::endl;
    std::cout << input["class"]["field"] << std::endl;
    std::cout << input["class"]["array"][0] << " " << input["class"]["array"][1] << " "<< input["class"]["array"][2] << std::endl;
    // One should usually force the type (c++ is strongly typed after all).
    input["number"].get<double>();
    input["class"]["field"].get<std::string>();

exit(0);
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

    std::ofstream file("he_data_lerp_1e4_splines.txt");
    file << std::setprecision(8) << std::scientific;
    // for (int i = 0; i < detectors.size(); i++)
    //     file << i << "\t";
    // file << "\n";

    double max = 0;
    
    for (int i = 0; i < H.size(); i++) {
        max = std::max(abs(spectrums[0][i]), max);
    }
    for (int i = 0; i < H.size(); i++) {
        spectrums[0][i] /= max;
    }
    for (int i = 0; i < H.size(); i++) {
        file << H[i] << "\t";
        for (int j = 0; j < detectors.size(); j++) {
            file << std::real(spectrums[j][i]) << "\t" << std::imag(spectrums[j][i]) << "\t";
        }
        file << "\n";
    }


    Profile::Print();

    return 0;
}