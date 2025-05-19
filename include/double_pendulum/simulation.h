#ifndef DOUBLE_PENDULUM_SIMULATION_H
#define DOUBLE_PENDULUM_SIMULATION_H

#include "double_pendulum/core.h"
#include <vector>
#include <cstddef> // For size_t

namespace double_pendulum {

// Runs the full simulation and returns all TwoPoints states
std::vector<TwoPoints> run_simulation_full(
    double L1, double L2, double M1, double M2,
    double p10, double p20, double o10, double o20,
    double T, double h);

// Runs a simulation and stores only the specified angle (phi1 or phi2) over time
void run_simulation_phi(
    double L1, double L2, double M1, double M2,
    double p10, double p20, double o10, double o20,
    double T, double h, bool use_phi1,
    std::vector<double>& phi_out); // Output vector for the angle
	
}

#endif // DOUBLE_PENDULUM_SIMULATION_H