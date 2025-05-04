#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <cstdlib>

<<<<<<< HEAD
// Platformfüggő beállítások
#if defined(_WIN32) || defined(_WIN64)
    #ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif

    #define popen _popen
    #define pclose _pclose
#endif

=======
>>>>>>> 20bb159f3871ee39fd1d2f26a7dd2334f4b1876f
const double g = 9.81;
const double T = 30.0;
const double h = 0.01;
const size_t N = static_cast<size_t>(T / h) + 1;
const int RES = 1000;
const double MIN_PHI = -3.1;
const double MAX_PHI = 3.1;

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

void simulate_phi(double p10, double p20, bool use_phi1, std::vector<double>& phi_out) {
    double L1 = 1.0, L2 = 1.0, M1 = 1.0, M2 = 1.0;
    double p1 = p10, p2 = p20, o1 = 0.0, o2 = 0.0;

    for (size_t i = 0; i < N; ++i) {
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

double find_first_flip_time(const std::vector<double>& phi) {
    for (size_t i = 1; i < N; ++i) {
        double prev = phi[i - 1];
        double curr = phi[i];
        double lower = std::min(prev, curr);
        double upper = std::max(prev, curr);

        int k_start = static_cast<int>(std::ceil((lower / M_PI - 1.0) / 2.0));
        int k_end = static_cast<int>(std::floor((upper / M_PI - 1.0) / 2.0));

        if (k_start <= k_end) {
            double target = (2 * k_start + 1) * M_PI;
            double ratio = (target - prev) / (curr - prev);
            return h * (i - 1 + ratio);
        }
    }
    return -1.0;
}
void plot_data(const char* title, const char* xlabel, const char* ylabel,
               const std::vector<std::vector<double>>& data,
               const std::string& filename_png) {
    FILE* gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        std::cerr << "Nem sikerült megnyitni a Gnuplotot.\n";
        return;
    }

    fprintf(gp,
        "set terminal pngcairo size 1000,800 enhanced font 'Helvetica,12'\n"
        "set output '%s'\n"
        "set title '%s' font \",14\"\n"
        "set xlabel '%s' font \",12\"\n"
        "set ylabel '%s' font \",12\"\n"
        "set xrange [%f:%f]\n"
        "set yrange [%f:%f]\n"
        "set pm3d map\n"
        "set size ratio -1\n"
        "set cbrange [-1:*]\n"
        "set palette defined ( \\\n"
        "    -1 '#ffffff', \\\n"
        "     0 '#440154', \\\n"
        "     1 '#414487', \\\n"
        "     2 '#2a788e', \\\n"
        "     3 '#22a884', \\\n"
        "     4 '#7ad151', \\\n"
        "     5 '#fde725')\n"
        "set key off\n"
        "splot '-' using 1:2:3 with image\n",
        filename_png.c_str(), title, xlabel, ylabel,
        MIN_PHI, MAX_PHI, MIN_PHI, MAX_PHI);

    for (int j = 0; j < RES; ++j) {
        double phi2 = MIN_PHI + (MAX_PHI - MIN_PHI) * j / (RES - 1);
        for (int i = 0; i < RES; ++i) {
            double phi1 = MIN_PHI + (MAX_PHI - MIN_PHI) * i / (RES - 1);
            fprintf(gp, "%f %f %f\n", phi1, phi2, data[i][j]);
        }
        fprintf(gp, "\n");  // az image plot elválasztónak szüksége van új sorra
    }

    fprintf(gp, "e\nunset output\n");
    pclose(gp);
}

int main() {
    std::vector<double> phi_arr(N);
    std::vector<std::vector<double>> flip1_data(RES, std::vector<double>(RES));
    std::vector<std::vector<double>> flip2_data(RES, std::vector<double>(RES));

    for (int i = 0; i < RES; ++i) {
        double phi1_0 = MIN_PHI + (MAX_PHI - MIN_PHI) * i / (RES - 1);
        for (int j = 0; j < RES; ++j) {
            double phi2_0 = MIN_PHI + (MAX_PHI - MIN_PHI) * j / (RES - 1);

            simulate_phi(phi1_0, phi2_0, true, phi_arr);
            flip1_data[i][j] = find_first_flip_time(phi_arr);

            simulate_phi(phi1_0, phi2_0, false, phi_arr);
            flip2_data[i][j] = find_first_flip_time(phi_arr);
        }
        std::cout << "Sor: " << i + 1 << " / " << RES << "\r" << std::flush;
    }

    std::cout << "\nKep mentes PNG-be...\n";

    plot_data("Elso atfordulas ideje: {/Symbol f}_1",
              "{/Symbol f}_1 (rad)", "{/Symbol f}_2 (rad)",
              flip1_data, "phi1.png");

    plot_data("Elso atfordulas ideje: {/Symbol f}_2",
              "{/Symbol f}_1 (rad)", "{/Symbol f}_2 (rad)",
              flip2_data, "phi2.png");

    return 0;
}
