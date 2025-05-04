#include <iostream>
#include <fstream>    
#include <cmath>
#include <string>
#include <vector>

const double g = 9.81; // nehézségi gyorsulás

struct TwoPoints {
    double time = 0.0;
    double x1, y1, v1;
    double x2, y2, v2;
    double phi1, phi2;
    double omega1, omega2;
};

TwoPoints converter(double L1, double L2, double p1, double p2, double o1, double o2) {
    TwoPoints coords;
    coords.x1 = L1 * sin(p1);
    coords.y1 = -L1 * cos(p1);
    coords.x2 = coords.x1 + L2 * sin(p2);
    coords.y2 = coords.y1 - L2 * cos(p2);
    coords.v1 = o1 * L1;
    coords.v2 = sqrt(pow(L1 * o1, 2) + pow(L2 * o2, 2) + 2 * L1 * L2 * o1 * o2 * cos(p1 - p2));
    coords.phi1 = p1;
    coords.phi2 = p2;
    coords.omega1 = o1;
    coords.omega2 = o2;
    return coords;
}

double do1dt(double L1, double L2, double M1, double M2, double p1, double p2, double o1, double o2) {
    return (-g * (2 * M1 + M2) * sin(p1) - M2 * g * sin(p1 - 2 * p2) -
            2 * sin(p1 - p2) * M2 * (pow(o2, 2) * L2 + pow(o1, 2) * L1 * cos(p1 - p2))) /
           (L1 * (2 * M1 + M2 - M2 * cos(2 * p1 - 2 * p2)));
}

double do2dt(double L1, double L2, double M1, double M2, double p1, double p2, double o1, double o2) {
    return (2 * sin(p1 - p2) * (pow(o1, 2) * L1 * (M1 + M2) + g * (M1 + M2) * cos(p1) +
                                pow(o2, 2) * L2 * M2 * cos(p1 - p2))) /
           (L2 * (2 * M1 + M2 - M2 * cos(2 * p1 - 2 * p2)));
}

std::vector<TwoPoints> simulation(double L1, double L2, double M1, double M2,
                                   double p10, double p20, double o10, double o20,
                                   double T, double h) {
    std::vector<TwoPoints> simulation_values;

    double p1 = p10, p2 = p20, o1 = o10, o2 = o20;
    int steps = static_cast<int>(T / h);

    for (int i = 0; i <= steps; ++i) {
        TwoPoints current = converter(L1, L2, p1, p2, o1, o2);
        current.time = i * h;
        simulation_values.push_back(current);

        // Runge–Kutta értékek
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

int main(int argc, char** argv) {
    if (argc != 11) {
        std::cerr << "You gave only " << argc - 1
                  << " parameters. Input should look like this: l1 l2 m1 m2 phi1_0 phi2_0 omega1_0 omega2_0 T h" << std::endl;
        return 1;
    }

    double l1 = std::stod(argv[1]), l2 = std::stod(argv[2]), m1 = std::stod(argv[3]), m2 = std::stod(argv[4]);
    double phi1_0 = std::stod(argv[5]), phi2_0 = std::stod(argv[6]);
    double omega1_0 = std::stod(argv[7]), omega2_0 = std::stod(argv[8]);
    double T = std::stod(argv[9]), h = std::stod(argv[10]);

    std::vector<TwoPoints> sim_values = simulation(l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h);

    std::ofstream output("double_pendulum_output.csv");

    if (output.is_open()) {
        output << "time,x1,y1,x2,y2,v1,v2,phi1,phi2,omega1,omega2\n";
        for (const auto& point : sim_values) {
            output << point.time << ","
                   << point.x1 << "," << point.y1 << ","
                   << point.x2 << "," << point.y2 << ","
                   << point.v1 << "," << point.v2 << ","
                   << point.phi1 << "," << point.phi2 << ","
                   << point.omega1 << "," << point.omega2 << "\n";
        }
<<<<<<< HEAD
        std::cout << "Results saved to file: double_pendulum_output.csv\n";
    } else {
        std::cerr << "Failed to open file for writing.\n";
=======
        std::cout << "Eredmények elmentve a double_pendulum_output.csv fájlba.\n";
    } else {
        std::cerr << "Nem sikerült megnyitni a fájlt írásra.\n";
>>>>>>> 20bb159f3871ee39fd1d2f26a7dd2334f4b1876f
    }

    return 0;
}
