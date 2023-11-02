#include "sin2pulse.h"
#include "maths.h"
#include "../maths/constants.h"

Sin2Pulse::Sin2Pulse(double cep, double w0, double N, double E0) 
    : CEP(cep), w0(w0), N(N), E0(E0)
{}
// vector potential
point3 Sin2Pulse::A(double t) const {
    // double duration = 2*pi*N / w0;
    // double phase = cep - w0*duration/2.;
    double phase = CEP - pi*N;
    point3 val{0, 0, 0};
    val.x = (E0/w0)*sin(0.5*w0*t/N)*sin(0.5*w0*t/N)*sin(w0*t + phase);
    return val;
}
// electric field

point3 Sin2Pulse::E(double t) const {
    // double duration = 2*pi*N / w0;
    // double phase = cep - w0*duration/2.;
    double phase = CEP - pi*N;
    point3 val{0, 0, 0};
    val.x = -E0*(sin(0.5*w0*t/N)*sin(0.5*w0*t/N)*cos(w0*t + phase)
              + cos(0.5*w0*t/N)*sin(0.5*w0*t/N)*sin(w0*t + phase) / N);
    return val;
}
// integral of vector potential
double Sin2Pulse::alpha(double t) const {
    // double duration = 2*pi*N / w0;
    // double phase = cep - w0*duration/2.;
    double phase = CEP - pi*N;
    return (E0/w0) * ( -2.*cos(phase)
                       -2.*(N*N-1.)*cos(w0*t + phase) + 
                        N*(N+1.)*cos((N-1.)*w0*t/N + phase) +
                        N*(N-1.)*cos((N+1.)*w0*t/N + phase))/(4.*(N*N - 1.)*w0);
}
// integral of vector potential square
double Sin2Pulse::beta(double t) const {
    double a = 0;
    double b = t;
    double duration = 2*pi*N / w0;
    double phase = CEP - w0*duration/2.;
    // double phase = CEP - pi*N;
    return E0*E0*(-6.*(2.*(a - b)*w0 
        - sin(2.*(w0*a + phase)) 
        + sin(2.*(w0*b + phase)))
        + N*(16.*(sin(w0*a/N) - sin(w0*b/N))
            + 2.*(sin(2.*b*w0/N) - sin(2.*a*w0/N))
            + 8.*(sin((2. + 1./N)*w0*b + 2.*phase) 
            - sin((2. + 1./N)*w0*a + 2.*phase)) / (2.*N + 1.)
            + (sin(2.*((1. - 1./N)*w0*a + phase))
            - sin(2.*((1. - 1./N)*w0*b + phase))) / (N - 1.)
            + (sin(2.*((1. + 1./N)*w0*a + phase))
            - sin(2.*((1. + 1./N)*w0*b + phase))) / (N + 1.)
            + 8.*(sin((2. - 1./N)*w0*b + 2.*phase) 
            - sin((2. - 1./N)*w0*a + 2.*phase)) / (2.*N - 1.)
        )) / (64.*w0*w0*w0);
}