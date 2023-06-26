#pragma once

#include <vector>
#include <complex>

namespace Spline {
    class Cubic {
        std::vector<std::complex<double>> a, b, c, d;
    public:
        void Initialize(const std::vector<double>& x, const std::vector<std::complex<double>>& y);
        std::complex<double> operator()(double x_s, const std::vector<double>& x) const;
    };
    void thomas_algorithm(  const std::vector<double>& a,
                            const std::vector<double>& b,
                            const std::vector<double>& c,
                            const std::vector<std::complex<double>>& d,
                            std::vector<std::complex<double>>& f);

};
