#pragma once

#include <vector>
#include <functional>
#include "../maths/constants.h"


//-------------------1d definite integrator using trapezoidal rule -------------/

template <typename T>
T Trapz(const std::vector<double>& xs, const std::vector<T>& f) {
    T sum = 0;
    for (int i = 1; i < xs.size(); i++) {
        sum += 0.5*(f[i] + f[i-1])*(xs[i] - xs[i-1]);
    }
    return sum;
}
template <typename T>
T Trapz(const std::vector<double>& xs, std::function<T(double)> f) {
    T sum = 0;
    for (int i = 1; i < xs.size(); i++) {
        sum += 0.5*(f(xs[i]) + f(xs[i-1]))*(xs[i] - xs[i-1]);
    }
    return sum;
}
template <typename T>
T Trapz(double xmin, double xmax, int nx, std::function<T(double)> f) {
    double dx = (xmax - xmin) / (nx - 1);
    T sum = 0;
    for (int i = 1; i < nx; i++) {
        double x = xmin + dx*i;
        sum += 0.5*(f(x) + f(x-dx))*dx;
    }
    return sum;
}

// ---------------2d definite intregator using trapezoidal rule -------------------------------/
template <typename T>
T Trapz2D(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<std::vector<T>>& f) {
    T sum = 0;
    for (int ix = 1; ix < xs.size(); ix++) {
        for (int iy = 1; iy < ys.size(); iy++) {
            double dx = xs[ix] - xs[ix - 1];
            double dy = ys[iy] - ys[iy - 1];
            sum += 0.25 * (f[ix][iy] + f[ix - 1][iy] + f[ix][iy - 1] + f[ix - 1][iy - 1]) * dx * dy;
        }
    }
    return sum;
}

template <typename T>
T Trapz2D(const std::vector<double>& xs, const std::vector<double>& ys, std::function<T(double, double)> f) {
    T sum = 0;
    for (int ix = 1; ix < xs.size(); ix++) {
        for (int iy = 1; iy < ys.size(); iy++) {
            double dx = xs[ix] - xs[ix - 1];
            double dy = ys[iy] - ys[iy - 1];
            sum += 0.25 * (f(xs[ix], ys[iy]) + f(xs[ix - 1], ys[iy]) + f(xs[ix], ys[iy - 1]) + f(xs[ix - 1], ys[iy - 1])) * dx * dy;
        }
    }
    return sum;
}

template <typename T>
T Trapz2D(double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::function<T(double, double)> f) {
    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);
    T sum = 0;
    for (int ix = 1; ix < nx; ix++) {
        for (int iy = 1; iy < ny; iy++) {
            double x = xmin + dx * ix;
            double y = ymin + dy * iy;
            sum += 0.25 * (f(x, y) + f(x - dx, y) + f(x, y - dy) + f(x - dx, y - dy)) * dx * dy;
        }
    }
    return sum;
}

//1d indefinite integrator using trapezoidal rule, sum is stored at each step
template <typename T>
std::vector<T> TrapzInd(const std::vector<double>& xs, const std::vector<T>& f) {
    std::vector<T> sum;
    sum.push_back(0);
    for (int i = 1; i < xs.size(); i++) {
        sum.push_back(sum[i - 1] + 0.5 * (f[i] + f[i - 1]) * (xs[i] - xs[i - 1]));
    }
    return sum;
}

template <typename T>
std::vector<T> TrapzInd(const std::vector<double>& xs, std::function<T(double)> f) {
    std::vector<T> sum;
    sum.push_back(0);
    for (int i = 1; i < xs.size(); i++) {
        sum.push_back(sum[i - 1] + 0.5 * (f(xs[i]) + f(xs[i - 1])) * (xs[i] - xs[i - 1]));
    }
    return sum;
}

template <typename T>
std::vector<T> TrapzInd(double xmin, double xmax, int nx, std::function<T(double)> f) {
    double dx = (xmax - xmin) / (nx - 1);
    std::vector<T> sum;
    sum.push_back(0);
    for (int i = 1; i < nx; i++) {
        double x = xmin + dx * i;
        sum.push_back(sum[i - 1] + 0.5 * (f(x) + f(x - dx)) * dx);
    }
    return sum;
}















