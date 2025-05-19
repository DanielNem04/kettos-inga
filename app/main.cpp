#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "double_pendulum/core.h"
#include "double_pendulum/simulation.h"
#include "double_pendulum/analysis.h"

using double_pendulum::TwoPoints;
using double_pendulum::run_simulation_full;
using double_pendulum::run_simulation_phi;
using double_pendulum::find_first_flip_time;

bool parse_args_csv(int argc, char** argv, double& l1, double& l2, double& m1, double& m2,
                    double& phi1_0, double& phi2_0, double& omega1_0, double& omega2_0,
                    double& T, double& h) {
    if (argc != 12) return false;
    try {
        l1 = std::stod(argv[2]); l2 = std::stod(argv[3]);
        m1 = std::stod(argv[4]); m2 = std::stod(argv[5]);
        phi1_0 = std::stod(argv[6]); phi2_0 = std::stod(argv[7]);
        omega1_0 = std::stod(argv[8]); omega2_0 = std::stod(argv[9]);
        T = std::stod(argv[10]); h = std::stod(argv[11]);
    } catch (...) { return false; }
    return true;
}

bool parse_args_flip_single(int argc, char** argv, double& l1, double& l2, double& m1, double& m2,
                            double& phi1_0, double& phi2_0, double& omega1_0, double& omega2_0,
                            double& T, double& h, int& use_phi1_int) {
    if (argc != 13) return false;
    try {
        l1 = std::stod(argv[2]); l2 = std::stod(argv[3]);
        m1 = std::stod(argv[4]); m2 = std::stod(argv[5]);
        phi1_0 = std::stod(argv[6]); phi2_0 = std::stod(argv[7]);
        omega1_0 = std::stod(argv[8]); omega2_0 = std::stod(argv[9]);
        T = std::stod(argv[10]); h = std::stod(argv[11]);
        use_phi1_int = std::stoi(argv[12]);
        if (use_phi1_int != 0 && use_phi1_int != 1) return false;
    } catch (...) { return false; }
    return true;
}

bool parse_args_flip_grid(int argc, char** argv, double& l1, double& l2, double& m1, double& m2,
                          double& omega1_0, double& omega2_0, double& T, double& h,
                          double& min_phi, double& max_phi, int& resolution) {
    if (argc != 13) return false;
    try {
        l1 = std::stod(argv[2]); l2 = std::stod(argv[3]);
        m1 = std::stod(argv[4]); m2 = std::stod(argv[5]);
        omega1_0 = std::stod(argv[6]); omega2_0 = std::stod(argv[7]);
        T = std::stod(argv[8]); h = std::stod(argv[9]);
        min_phi = std::stod(argv[10]); max_phi = std::stod(argv[11]);
        resolution = std::stoi(argv[12]);
        if (resolution <= 0 || min_phi >= max_phi) return false;
    } catch (...) { return false; }
    return true;
}

std::string construct_filename(const std::string& prefix, double l1, double l2, double m1, double m2,
                               double omega1_0, double omega2_0, double T, double h,
                               double min_phi, double max_phi, int res) {
    std::ostringstream oss;
    oss << prefix
        << "_l1_" << l1 << "_l2_" << l2
        << "_m1_" << m1 << "_m2_" << m2
        << "_o1_" << omega1_0 << "_o2_" << omega2_0
        << "_T_" << T << "_h_" << h
        << "_phi_" << min_phi << "_to_" << max_phi
        << "_res_" << res << ".csv";
    return oss.str();
}

void write_csv(const std::string& filename, const std::vector<std::vector<double>>& data, double min_phi, double max_phi) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    int size = static_cast<int>(data.size());
    double step = (max_phi - min_phi) / (size > 1 ? size - 1 : 1.0);

    file << "phi1,phi2,flip_time\n";
    for (int i = 0; i < size; ++i) {
        double phi1 = min_phi + step * i;
        for (int j = 0; j < size; ++j) {
            double phi2 = min_phi + step * j;
            file << phi1 << "," << phi2 << "," << data[i][j] << "\n";
        }
    }
    std::cout << "Saved CSV: " << filename << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <mode> [parameters...]\n";
        std::cerr << "Modes:\n"
                  << "  csv l1 l2 m1 m2 phi1_0 phi2_0 omega1_0 omega2_0 T h\n"
                  << "  flip_single l1 l2 m1 m2 phi1_0 phi2_0 omega1_0 omega2_0 T h use_phi1(0 or 1)\n"
                  << "  flip_grid l1 l2 m1 m2 omega1_0 omega2_0 T h min_phi max_phi resolution\n";
        return 1;
    }

    std::string mode = argv[1];

    if (mode == "csv") {
        double l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h;
        if (!parse_args_csv(argc, argv, l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h)) {
            std::cerr << "Invalid CSV mode arguments.\n";
            return 1;
        }
        auto sim_values = run_simulation_full(l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h);
        std::ofstream out("../../app/simulation_output.csv");
        out << "time,x1,y1,x2,y2,v1,v2,phi1,phi2,omega1,omega2\n";
        for (const auto& pt : sim_values) {
            out << pt.time << "," << pt.x1 << "," << pt.y1 << "," << pt.x2 << "," << pt.y2 << ","
                << pt.v1 << "," << pt.v2 << "," << pt.phi1 << "," << pt.phi2 << ","
                << pt.omega1 << "," << pt.omega2 << "\n";
        }
        std::cout << "Saved simulation_output.csv\n";
    }
    else if (mode == "flip_single") {
        double l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h;
        int use_phi1_int;
        if (!parse_args_flip_single(argc, argv, l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h, use_phi1_int)) {
            std::cerr << "Invalid flip_single arguments.\n";
            return 1;
        }
        std::vector<double> phi_values;
        run_simulation_phi(l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h, use_phi1_int == 1, phi_values);
        double flip_time = find_first_flip_time(phi_values, h);
        std::cout << "First flip time (" << (use_phi1_int ? "phi1" : "phi2") << "): " << flip_time << std::endl;
    }
    else if (mode == "flip_grid") {
        double l1, l2, m1, m2, omega1_0, omega2_0, T, h, min_phi, max_phi;
        int resolution;
        if (!parse_args_flip_grid(argc, argv, l1, l2, m1, m2, omega1_0, omega2_0, T, h, min_phi, max_phi, resolution)) {
            std::cerr << "Invalid flip_grid arguments.\n";
            return 1;
        }

        std::vector<std::vector<double>> flip1_data(resolution, std::vector<double>(resolution));
        std::vector<std::vector<double>> flip2_data(resolution, std::vector<double>(resolution));
        std::vector<double> phi_arr;
        double step = (max_phi - min_phi) / (resolution > 1 ? resolution - 1 : 1);

        for (int i = 0; i < resolution; ++i) {
            double phi1_0 = min_phi + step * i;
            for (int j = 0; j < resolution; ++j) {
                double phi2_0 = min_phi + step * j;
                phi_arr.clear();
                run_simulation_phi(l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h, true, phi_arr);
                flip1_data[i][j] = find_first_flip_time(phi_arr, h);
                phi_arr.clear();
                run_simulation_phi(l1, l2, m1, m2, phi1_0, phi2_0, omega1_0, omega2_0, T, h, false, phi_arr);
                flip2_data[i][j] = find_first_flip_time(phi_arr, h);
            }
            std::cout << "Progress: " << (i + 1) * 100 / resolution << "%\r" << std::flush;
        }
        std::cout << std::endl;

        std::string file1 = "../../app/" + construct_filename("phi1_flip", l1, l2, m1, m2, omega1_0, omega2_0, T, h, min_phi, max_phi, resolution);
        std::string file2 = "../../app/" + construct_filename("phi2_flip", l1, l2, m1, m2, omega1_0, omega2_0, T, h, min_phi, max_phi, resolution);

        write_csv(file1, flip1_data, min_phi, max_phi);
        write_csv(file2, flip2_data, min_phi, max_phi);
    }
    else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        return 1;
    }
    return 0;
}