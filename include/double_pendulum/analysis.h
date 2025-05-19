#ifndef DOUBLE_PENDULUM_ANALYSIS_H
#define DOUBLE_PENDULUM_ANALYSIS_H

#include <vector>
#include <cmath> // For std::ceil, std::floor, std::abs
#include <cstddef> // For size_t

// Ensure M_PI is defined, if not already present
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace double_pendulum {

// Finds the time of the first "flip" (crossing (2k+1)*pi)
// Returns the time of the flip, or -1.0 if no flip is found within the data.
double find_first_flip_time(const std::vector<double>& phi, double h);

}


#endif // DOUBLE_PENDULUM_ANALYSIS_H