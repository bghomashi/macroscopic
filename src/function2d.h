#pragma once

#include "maths.h"

#include <vector>

struct Function2D {
    std::vector<complex> range;
    std::vector<double> xs;
    std::vector<double> ys;

    Function2D(){}
    Function2D(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<complex>& f) : xs(xs), ys(ys) {
        range = f;
    }
    Function2D(const std::vector<double>& xs, const std::vector<double>& ys) : xs(xs), ys(ys), range(xs.size()*ys.size()) {}

    inline int findNearestLowerX(double x) const {
        for (int i = 0; i < xs.size()-1; i++)
            if (x >= xs[i] && x < xs[i+1])
                return i;
        return xs.size() - 1;
    }
    inline int findNearestLowerY(double y) const  {
        for (int i = 0; i < ys.size()-1; i++)
            if (y >= ys[i] && y < ys[i+1])
                return i;
        return ys.size() - 1;
    }
    inline complex operator ()(double x, double y) const {
        int ix = findNearestLowerX(x);
        int iy = findNearestLowerY(y);
        
        return get(ix, iy);
    }
    inline complex& operator ()(double x, double y) {
        int ix = findNearestLowerX(x);
        int iy = findNearestLowerY(y);
        
        return get(ix, iy);
    }
    inline complex get(int ix, int iy) const {
        int index = ix + iy*xs.size();
        return range[index];
    }
    inline complex& get(int ix, int iy) {
        int index = ix + iy*xs.size();
        return range[index];
    }

    static void Fourier(const Function2D& f_in, Function2D& f_out) {
        for (int ikx = 0; ikx < f_out.xs.size(); ikx++) {
            for (int iky = 0; iky < f_out.ys.size(); iky++) {
                double kx = f_out.xs[ikx];
                double ky = f_out.ys[iky];

                auto& Fxy = f_out.get(ikx, iky);
                Fxy = 0;
                for (int ix = 1; ix < f_in.xs.size(); ix++) {
                    double dx = f_in.xs[ix] - f_in.xs[ix-1];
                    double x0 = f_in.xs[ix];
                    double x1 = f_in.xs[ix-1];
                    complex I0y = 0;
                    complex I1y = 0;

                    for (int iy = 1; iy < f_in.ys.size(); iy++) {
                        double dy = f_in.ys[iy] - f_in.ys[iy-1];
                        double y0 = f_in.ys[iy];
                        double y1 = f_in.ys[iy-1];

                        auto F00 = f_in.get(ix  , iy  ) * exp(-I*(kx*x0 + ky*y0));
                        auto F10 = f_in.get(ix-1, iy  ) * exp(-I*(kx*x1 + ky*y0));
                        auto F01 = f_in.get(ix  , iy-1) * exp(-I*(kx*x0 + ky*y1));
                        auto F11 = f_in.get(ix-1, iy-1) * exp(-I*(kx*x1 + ky*y1));

                        I0y += 0.5 * (F00 + F01) * dy;
                        I1y += 0.5 * (F10 + F11) * dy;
                    }
                    Fxy += 0.5 * (I0y + I1y) * dx;
                }
            }
        }
    }
};