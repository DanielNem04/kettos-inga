#include "double_pendulum/analysis.h"
#include <algorithm> // For std::min, std::max

namespace double_pendulum {


double find_first_flip_time(const std::vector<double>& phi, double h) {
    for (size_t i = 1; i < phi.size(); ++i) {
        double prev = phi[i - 1];
        double curr = phi[i];

        // We need to find integer k such that (2k+1)*pi is between prev and curr vals (or curr and prev).
        // This is equivalent to checking if prev and curr are on opposite sides of any odd multiple of pi.

        double lower_bound = std::min(prev, curr);
        double upper_bound = std::max(prev, curr);

        int k_start = static_cast<int>(std::ceil((lower_bound / M_PI - 1.0) / 2.0));
        int k_end = static_cast<int>(std::floor((upper_bound / M_PI - 1.0) / 2.0));


        for (int k = k_start; k <= k_end; ++k) {
            double target = (2 * k + 1) * M_PI;

             // Check if the target is strictly between prev and curr.
            if ((prev < target && target < curr) || (curr < target && target < prev)) {
                // Linear interpolation to find the time of crossing not necessary in our case, but why not?

                double ratio = (target - prev) / (curr - prev);
                return h * (static_cast<double>(i - 1) + ratio);
            }
        }
    }
    return -1.0; // No flip found (i thought it is convinient for plotting )
}

}