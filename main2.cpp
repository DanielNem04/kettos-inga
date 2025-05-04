#include <cmath>
#include <cstddef>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


extern "C" {

const double g = 9.81;

// RK4 deriv
double do1dt(double L1, double L2, double M1, double M2,
             double p1, double p2, double o1, double o2) {
    return (-g * (2 * M1 + M2) * sin(p1) - M2 * g * sin(p1 - 2 * p2) -
            2 * sin(p1 - p2) * M2 * (pow(o2, 2) * L2 + pow(o1, 2) * L1 * cos(p1 - p2))) /
           (L1 * (2 * M1 + M2 - M2 * cos(2 * p1 - 2 * p2)));
}

double do2dt(double L1, double L2, double M1, double M2,
             double p1, double p2, double o1, double o2) {
    return (2 * sin(p1 - p2) * (pow(o1, 2) * L1 * (M1 + M2) + g * (M1 + M2) * cos(p1) +
                                pow(o2, 2) * L2 * M2 * cos(p1 - p2))) /
           (L2 * (2 * M1 + M2 - M2 * cos(2 * p1 - 2 * p2)));
}

// RK4 szimu
void simulate_phi(double L1, double L2, double M1, double M2,
                  double p10, double p20, double o10, double o20,
                  double T, double h, int use_phi1,
                  double* phi_out, double* time_out, size_t length) {
    double p1 = p10, p2 = p20, o1 = o10, o2 = o20;

    for (size_t i = 0; i < length; ++i) {
        time_out[i] = i * h;
        phi_out[i] = use_phi1 ? p1 : p2;

        double k1_o1 = do1dt(L1, L2, M1, M2, p1, p2, o1, o2);
        double k1_o2 = do2dt(L1, L2, M1, M2, p1, p2, o1, o2);

        double k2_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2,
                             o1 + h * k1_o1 / 2, o2 + h * k1_o2 / 2);
        double k2_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2,
                             o1 + h * k1_o1 / 2, o2 + h * k1_o2 / 2);

        double k3_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2,
                             o1 + h * k2_o1 / 2, o2 + h * k2_o2 / 2);
        double k3_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2,
                             o1 + h * k2_o1 / 2, o2 + h * k2_o2 / 2);

        double k4_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1, p2 + h * o2,
                             o1 + h * k3_o1, o2 + h * k3_o2);
        double k4_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1, p2 + h * o2,
                             o1 + h * k3_o1, o2 + h * k3_o2);

        o1 += h * (k1_o1 + 2 * k2_o1 + 2 * k3_o1 + k4_o1) / 6.0;
        o2 += h * (k1_o2 + 2 * k2_o2 + 2 * k3_o2 + k4_o2) / 6.0;

        p1 += h * o1;
        p2 += h * o2;
    }
}

// Átfordulás keresése: első (2k+1)*pi-n való áthaladás
double find_first_flip_time(const double* time, const double* phi, size_t length) {
    for (size_t i = 1; i < length; ++i) {
        double prev = phi[i - 1];
        double curr = phi[i];
        double lower = (prev < curr) ? prev : curr;
        double upper = (prev > curr) ? prev : curr;

        int k_start = static_cast<int>(std::ceil((lower / M_PI - 1.0) / 2.0));
        int k_end = static_cast<int>(std::floor((upper / M_PI - 1.0) / 2.0));

        if (k_start <= k_end) {
            double target = (2 * k_start + 1) * M_PI;
            double ratio = (target - prev) / (curr - prev);
            return time[i - 1] + ratio * (time[i] - time[i - 1]);
        }
    }
    return -1.0;
}

} // extern "C"
