#include "double_pendulum/core.h"



namespace double_pendulum {

const double g = 9.81; // g const

TwoPoints converter(double L1, double L2, double p1, double p2, double o1, double o2) {
    TwoPoints coords;
    coords.x1 = L1 * std::sin(p1);
    coords.y1 = -L1 * std::cos(p1);
    coords.x2 = coords.x1 + L2 * std::sin(p2);
    coords.y2 = coords.y1 - L2 * std::cos(p2);
    coords.v1 = o1 * L1;
    coords.v2 = std::sqrt(std::pow(L1 * o1, 2) + std::pow(L2 * o2, 2) + 2 * L1 * L2 * o1 * o2 * std::cos(p1 - p2));
    coords.phi1 = p1;
    coords.phi2 = p2;
    coords.omega1 = o1;
    coords.omega2 = o2;
    return coords;
}

double do1dt(double L1, double L2, double M1, double M2, double p1, double p2, double o1, double o2) {
    return (-g * (2 * M1 + M2) * std::sin(p1) - M2 * g * std::sin(p1 - 2 * p2) -
            2 * std::sin(p1 - p2) * M2 * (std::pow(o2, 2) * L2 + std::pow(o1, 2) * L1 * std::cos(p1 - p2))) /
           (L1 * (2 * M1 + M2 - M2 * std::cos(2 * p1 - 2 * p2)));
}

double do2dt(double L1, double L2, double M1, double M2, double p1, double p2, double o1, double o2) {
    return (2 * std::sin(p1 - p2) * (std::pow(o1, 2) * L1 * (M1 + M2) + g * (M1 + M2) * std::cos(p1) +
                                std::pow(o2, 2) * L2 * M2 * std::cos(p1 - p2))) /
           (L2 * (2 * M1 + M2 - M2 * std::cos(2 * p1 - 2 * p2)));
}

}