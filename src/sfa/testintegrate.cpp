// test_integration.cpp

#include "integrate.h"
#include <iostream>
#include <fstream>
#include <cmath>

int main() {
    // Test 1D definite integrator using trapezoidal rule
    std::vector<double> xs = {0.0, 1.0, 2.0, 3.0, 3.5};
    std::vector<double> ys = {0.0, 1.0, 2.0, 3.0, 3.5};
    
    std::function<double(double)> f1d = [](double x) { return x * x; };
    std::function<double(double,double)> f2d = [](double x, double y) { return x * x + y * y; };

    auto result1D = Trapz(xs, f1d);
    auto result2D = Trapz2D(xs, ys, f2d);

    std::cout << "1D Definite Integral: " << result1D << std::endl;
    std::cout << "2D Definite Integral: " << result2D << std::endl;


    // Test 1D indefinite integrator using trapezoidal rule
    int nx = 100;
    std::function<double(double)> f1dIndefinite = [](double x) { return sin(x); };

    std::vector<double> xsIndefinite(nx);
    for (int i = 0; i < nx; ++i) {
        xsIndefinite[i] = i * 0.1;
    }

    std::vector<double> sum1DIndefinite = TrapzInd(xsIndefinite, f1dIndefinite);

    std::ofstream outFile1DIndefinite("result_1d_indefinite.txt");
    if (outFile1DIndefinite.is_open()) {
        for (size_t i = 0; i < sum1DIndefinite.size(); ++i) {
            outFile1DIndefinite << xsIndefinite[i] << " " << sum1DIndefinite[i] << std::endl;
        }
        outFile1DIndefinite.close();
    }

    return 0;
}