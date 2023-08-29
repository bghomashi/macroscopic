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
std::vector<std::pair<double, double>> detectors;
std::vector<cvector> spectrums(detectors.size());

int main() {
    // HOW TO READ JSON
    nlohmann::json input;                           // create json-object
    std::ifstream in_file_stream("input.json");     // open the (text) file
    in_file_stream >> input;                        // read the file into (>>) the json-object


    Macroscopic macroscopic;
    
    
    points _points;
    //set up points
    if (input["samples"]["custom"] == "" && input["samples"]["preset"] == "") {
        std::cout << "No sample specified" << std::endl;
        exit(1);
    }
    if (input["samples"]["custom"] != "") {
        //still need to impliment how to read in custom points
    }
    else if (input["sample"]["preset"] == "cylindrical_gas_jet")
    {
        _points.CylindricalGasJet(
            input["samples"]["cylindrical_gas_jet"]["density_cm3"],
            input["samples"]["cylindrical_gas_jet"]["sigma_um"],
            input["samples"]["cylindrical_gas_jet"]["radius_um"],
            input["samples"]["cylindrical_gas_jet"]["length_um"],
            input["samples"]["cylindrical_gas_jet"]["number"]
        );
        _points.cutoff = input["samples"]["cylindrical_gas_jet"]["cutoff"];
    }
    else if (input["samples"]["preset"] == "line")
    {
        //need to impliment
    }
    else if (input["samples"]["preset"] == "square_lattice")
    {
        //need to impliment
    }
    else if (input["samples"]["preset"] == "cube_lattice")
    {
        //need to impliment
    }
    else 
    {
        std::cout << "Invalid sample preset" << std::endl;
        exit(1);
    }
    double x_shift = input["samples"]["shift"]["x_um"];
    double y_shift = input["samples"]["shift"]["y_um"];
    double z_shift = input["samples"]["shift"]["z_um"];
    //shift all the points  by the specified amount
    for (int i = 0; i < _points._cells.size(); i++) {
        _points._cells[i].pos.x += x_shift;
        _points._cells[i].pos.y += y_shift;
        _points._cells[i].pos.z += z_shift;
    }


    //set up laser
    Laser _laser(
        input["laser"]["wavelength_nm"],
        input["laser"]["beam_waist_um"],
        input["laser"]["peak_I0_wcm2"],
        input["laser"]["porras"]
    );
    

    //set up spectrum
    filename = input["spectrum"]["filename"];

    //set up detectors
    if(input["detector"]["preset"] == "on_axis")
    {
        detectors.push_back({0, 0});
    }
    else if(input["detector"]["preset"] == "point")
    {
        detectors.push_back({input["detector"]["point"]["theta_pi"], input["detector"]["point"]["phi_pi"]});
    }
    else if(input["detector"]["preset"] == "custom")
    {
        //need to impliment
    }
    else if(input["detector"]["preset"] == "circular")
    {
        //need to impliment
    }
    else if(input["detector"]["preset"] == "spherical")
    {
        //need to impliment
    }
    else if(input["detector"]["preset"] == "linear")
    {
        //need to impliment
    }
    else if(input["detector"]["preset"] == "cartesian")
    {
        //need to impliment
    }
    else
    {
        std::cout << "Invalid detector preset" << std::endl;
        exit(1);
    } 

    //outputfile 
    std::string outputfile = input["output"];
    std::ofstream file(outputfile);
    file << std::setprecision(8) << std::scientific;
    file << "H\t" <<  "theta\t" << "phi\t" << "real\t" << "imag\t" << "\n";

    macroscopic.Initialize(
        _laser,
        _points,
        filename
    );
    auto H = macroscopic.Frequencies();
    for (int i = 0; i < detectors.size(); i++) {
        auto& d = detectors[i];
        double td = d.first;
        double pd = d.second;

        Profile::Push("Spectrum");
        //output spectrum to out out file
        auto spectrum = macroscopic.Spectrum(td, pd);
        for (int i = 0; i < H.size(); i++) {
            file << H[i] << "\t" << td << "\t" << pd << "\t" << std::real(spectrum[i]) << "\t" << std::imag(spectrum[i]) << "\n";
        }
        Profile::Pop("Spectrum");
    }





    Profile::Print();

    return 0;
}