#pragma once

#include "sfa.h"
#include "../maths/constants.h"

struct Sin2Pulse : Pulse {
    double CEP;
    double N;
    double w0;
    double E0;

    Sin2Pulse(double cep, double w0, double N, double E0);
    point3 A(double t) const;
    point3 E(double t) const;
    double alpha(double t) const;
    double beta(double t) const;
};