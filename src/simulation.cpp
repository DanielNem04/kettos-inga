#include "double_pendulum/simulation.h"
#include "double_pendulum/core.h" // Include core.h for do1dt, do2dt, converter

namespace double_pendulum {


std::vector<TwoPoints> run_simulation_full(
    double L1, double L2, double M1, double M2,
    double p10, double p20, double o10, double o20,
    double T, double h) {

    std::vector<TwoPoints> simulation_values;

    double p1 = p10, p2 = p20, o1 = o10, o2 = o20;
    int steps = static_cast<int>(T / h);
    if (steps < 0) steps = 0; //  potential negative T

    simulation_values.reserve(steps + 1); // this was its more optimized

    for (int i = 0; i <= steps; ++i) {
        TwoPoints current = converter(L1, L2, p1, p2, o1, o2);
        current.time = i * h;
        simulation_values.push_back(current);

        // Runge–Kutta 4th order
        double k1_o1 = do1dt(L1, L2, M1, M2, p1, p2, o1, o2);
        double k1_o2 = do2dt(L1, L2, M1, M2, p1, p2, o1, o2);

        double k2_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k1_o1 / 2, o2 + h * k1_o2 / 2);
        double k2_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k1_o1 / 2, o2 + h * k1_o2 / 2);

        double k3_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k2_o1 / 2, o2 + h * k2_o2 / 2);
        double k3_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k2_o1 / 2, o2 + h * k2_o2 / 2);

        double k4_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1, p2 + h * o2, o1 + h * k3_o1, o2 + h * k3_o2);
        double k4_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1, p2 + h * o2, o1 + h * k3_o1, o2 + h * k3_o2);

        o1 += h * (k1_o1 + 2 * k2_o1 + 2 * k3_o1 + k4_o1) / 6.0;
        o2 += h * (k1_o2 + 2 * k2_o2 + 2 * k3_o2 + k4_o2) / 6.0;

        p1 += h * o1;
        p2 += h * o2;
    }

    return simulation_values;
}
// same but just phi vals are returned
void run_simulation_phi(
    double L1, double L2, double M1, double M2,
    double p10, double p20, double o10, double o20,
    double T, double h, bool use_phi1,
    std::vector<double>& phi_out) {

    double p1 = p10, p2 = p20, o1 = o10, o2 = o20;
    int steps = static_cast<int>(T / h);
    if (steps < 0) steps = 0; 

    phi_out.resize(steps + 1); 

    for (int i = 0; i <= steps; ++i) {
        phi_out[i] = use_phi1 ? p1 : p2;

        double k1_o1 = do1dt(L1, L2, M1, M2, p1, p2, o1, o2);
        double k1_o2 = do2dt(L1, L2, M1, M2, p1, p2, o1, o2);

        double k2_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k1_o1 / 2, o2 + h * k1_o2 / 2);
        double k2_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k1_o1 / 2, o2 + h * k1_o2 / 2);

        double k3_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k2_o1 / 2, o2 + h * k2_o2 / 2);
        double k3_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1 / 2, p2 + h * o2 / 2, o1 + h * k2_o1 / 2, o2 + h * k2_o2 / 2);

        double k4_o1 = do1dt(L1, L2, M1, M2, p1 + h * o1, p2 + h * o2, o1 + h * k3_o1, o2 + h * k3_o2);
        double k4_o2 = do2dt(L1, L2, M1, M2, p1 + h * o1, p2 + h * o2, o1 + h * k3_o1, o2 + h * k3_o2);

        o1 += h * (k1_o1 + 2 * k2_o1 + 2 * k3_o1 + k4_o1) / 6.0;
        o2 += h * (k1_o2 + 2 * k2_o2 + 2 * k3_o2 + k4_o2) / 6.0;

        p1 += h * o1;
        p2 += h * o2;
    }
}


}