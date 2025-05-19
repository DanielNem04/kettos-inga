#ifndef DOUBLE_PENDULUM_CORE_H
#define DOUBLE_PENDULUM_CORE_H

#include <cmath> // For sin, cos, pow, sqrt
#include <vector>

// Ensure M_PI is defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace double_pendulum {


extern const double g; // g constan

struct TwoPoints {
    double time = 0.0;
    double x1, y1, v1;
    double x2, y2, v2;
    double phi1, phi2;
    double omega1, omega2;
};

// Converts angular to XYZ
TwoPoints converter(double L1, double L2, double p1, double p2, double o1, double o2);

// Calculates the derivative of omega1
double do1dt(double L1, double L2, double M1, double M2, double p1, double p2, double o1, double o2);

// Calculates the derivative of omega2
double do2dt(double L1, double L2, double M1, double M2, double p1, double p2, double o1, double o2);



}
#endif // DOUBLE_PENDULUM_CORE_H