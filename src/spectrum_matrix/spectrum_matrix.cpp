#include "spectrum_matrix.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "utility/profiler.h"

std::vector<std::complex<double>> SpectrumMatrix::Lerp(double intensity) const {
    std::vector<std::complex<double>> out_spectrum(frequencies.size());

    // find bounding intensities
    int intensity_index;
    for (intensity_index = 0; intensity_index < intensities.size()-1; intensity_index++) {
        if (intensity <= intensities[intensity_index] && intensity > intensities[intensity_index + 1])
            break;
    }

    for (int i = 0; i < frequencies.size(); i++) {
        std::complex<double> a = get(intensity_index+1, i);     // smaller one
        std::complex<double> b = get(intensity_index, i);       // larger

        // out_spectrum[i] = b;
        out_spectrum[i] = lerp( intensities[intensity_index + 1], intensities[intensity_index],
                                a, b,
                                intensity);

        // out_spectrum[i] = splines[i](intensity, intensities);
    }

    return out_spectrum;
}

SpectrumMatrix SpectrumMatrix::Load(const std::string& filename) {
    SpectrumMatrix matrix;
    std::string line;
    std::ifstream file(filename);
    if (!file.is_open())
        return matrix;

    // ---------------- first line in file contains frequencies
    std::getline(file, line);
    std::stringstream ss(line);
    double f;
    while (ss >> f) {
        matrix.frequencies.push_back(f);
    }
    // ---------------------------------------------------------
    // ---------------- now the first number contains the
    // ---------------- intensity and the rest is the HHG (complex)
    Profile::Push("load data");
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        matrix.intensities.push_back(0);
        ss >> matrix.intensities.back();
        double real, imag;
        while (ss >> real >> imag) {
            matrix.data.push_back({ real, imag });
        }
    }
    Profile::Pop("load data");
    for (auto& i : matrix.intensities)
        i *= Wcm2ToAU;
    // ---------------------------------------------------------
    Profile::Push("generate splines");
    size_t num_freq = matrix.frequencies.size();
    size_t num_inten = matrix.intensities.size();
    matrix.splines.resize(num_freq);
    cvector temp(num_inten);
    for (int i = 0; i < num_freq; i++) {
        for (int j = 0; j < num_inten; j++)
            temp[j] = matrix.get(j,i);

        // matrix.splines[i].Initialize(matrix.intensities, temp);
    }
    Profile::Pop("generate splines");
    Profile::Print();
    return matrix;
}





const std::vector<double>& SpectrumMatrix::Frequencies() const {
    return frequencies;
}
const std::vector<double>& SpectrumMatrix::Intensities() const {
    return intensities;
}
std::complex<double>& SpectrumMatrix::get(int i, int f) {
    return data[f + i * frequencies.size()];
}
std::complex<double> SpectrumMatrix::get(int i, int f) const {
    return data[f + i * frequencies.size()];
}
std::complex<double> SpectrumMatrix::lerp(double x0, double x1, std::complex<double> y0, std::complex<double> y1, double x) const {
    return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}