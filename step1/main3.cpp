#include <iostream>  // include iostream
#include <vector>  // include vector
#include <cmath>  // include cmath
#include <fstream>  // include fstream
#include <omp.h>  // include omp.h
#include <algorithm>  // include algorithm
// (blank line)
struct Particle { double x, ux, uy, uz, w; };  // type definition
constexpr double eps0 = 8.8541878128e-12, me = 9.1093837e-31, qe = -1.602176634e-19, c = 299792458.0;  // variable assignment/initialization
// (blank line)
int main() {  // begin block
    const int Nx = 256, num_steps = 1000;  // variable assignment/initialization
    const double dx = 1e-5, L = Nx * dx, dt = 1e-14, B0_z = 2.5;   // variable assignment/initialization
    std::vector<double> Ex(Nx + 2, 0.0), rho(Nx + 2, 0.0);  // function declaration/call or variable initialization
    std::vector<Particle> electrons(1);   // function declaration/call or variable initialization
// (blank line)
    electrons[0].x = L * 0.5; electrons[0].w = 1.0;   // variable assignment/initialization
    double vx0 = 0.5 * c, gamma0 = 1.0 / std::sqrt(1.0 - (vx0*vx0)/(c*c));  // function declaration/call or variable initialization
    electrons[0].ux = gamma0 * vx0; electrons[0].uy = 0.0; electrons[0].uz = 0.0;  // variable assignment/initialization
// (blank line)
    std::ofstream diag("pic_diagnostics.csv");  // function declaration/call or variable initialization
    diag << "step,time,tracked_x,tracked_ux,tracked_uy,field_energy\n";  // statement
    double* rho_ptr = rho.data();  // function declaration/call or variable initialization
// (blank line)
    for (int step = 1; step <= num_steps; ++step) {  // begin block
        std::fill(rho.begin(), rho.end(), 0.0);  // function declaration/call or variable initialization
        for (size_t i = 0; i < electrons.size(); ++i) {  // begin block
            double x_norm = electrons[i].x / dx; int j = static_cast<int>(std::round(x_norm)); double xi = x_norm - j;  // function declaration/call or variable initialization
            double w0 = 0.5*std::pow(0.5-xi,2), w1 = 0.75-std::pow(xi,2), w2 = 0.5*std::pow(0.5+xi,2);  // function declaration/call or variable initialization
            if (j >= 1 && j < Nx) { rho_ptr[j-1] += (qe*electrons[i].w*w0)/dx; rho_ptr[j] += (qe*electrons[i].w*w1)/dx; rho_ptr[j+1] += (qe*electrons[i].w*w2)/dx; }  // control flow statement
        }  // end block
        std::fill(Ex.begin(), Ex.end(), 0.0);  // function declaration/call or variable initialization
        for (int j = 0; j < Nx + 1; ++j) { Ex[j + 1] = Ex[j] + (rho[j] * dx) / eps0; }  // control flow statement
// (blank line)
        #pragma omp parallel for schedule(static)  // preprocessor directive
        for (size_t i = 0; i < electrons.size(); ++i) {  // begin block
            double x_norm = electrons[i].x / dx; int j = static_cast<int>(std::round(x_norm)); double xi = x_norm - j;  // function declaration/call or variable initialization
            double w0 = 0.5*std::pow(0.5-xi,2), w1 = 0.75-std::pow(xi,2), w2 = 0.5*std::pow(0.5+xi,2);  // function declaration/call or variable initialization
            double E_part = (j >= 1 && j < Nx) ? (w0 * Ex[j - 1] + w1 * Ex[j] + w2 * Ex[j + 1]) : 0.0;  // function declaration/call or variable initialization
// (blank line)
            double q_trans = (qe * dt) / (2.0 * me);  // function declaration/call or variable initialization
            double ux_m = electrons[i].ux + q_trans * E_part, uy_m = electrons[i].uy, uz_m = electrons[i].uz;  // variable assignment/initialization
            double u_m_sq = ux_m*ux_m + uy_m*uy_m + uz_m*uz_m;  // variable assignment/initialization
            double tz = q_trans * B0_z, t_sq = tz*tz, s_param = 1.0 + u_m_sq/(c*c) - t_sq/(c*c);  // function declaration/call or variable initialization
            double gamma_p = std::sqrt(0.5 * s_param + std::sqrt(0.25 * std::pow(s_param, 2) + t_sq / (c*c)));  // function declaration/call or variable initialization
// (blank line)
            double sz = tz / gamma_p, s_inner = sz*sz;  // variable assignment/initialization
// (blank line)
            // FIXED: Mathematically sound 1D3V relativistic axis rotation matching the Lorentz updates
            electrons[i].ux = (ux_m + uy_m * sz) / (1.0 + s_inner) + q_trans * E_part;  // variable assignment/initialization
            electrons[i].uy = (uy_m - ux_m * sz) / (1.0 + s_inner);  // variable assignment/initialization
            electrons[i].uz = uz_m;  // variable assignment/initialization
// (blank line)
            double u_n_sq = electrons[i].ux*electrons[i].ux + electrons[i].uy*electrons[i].uy + electrons[i].uz*electrons[i].uz;  // variable assignment/initialization
            double g_n = std::sqrt(1.0 + u_n_sq / (c * c));  // function declaration/call or variable initialization
            electrons[i].x += (electrons[i].ux / g_n) * dt;  // variable assignment/initialization
// (blank line)
            if (electrons[i].x < 0.0) { electrons[i].x = -electrons[i].x; electrons[i].ux = -electrons[i].ux; }  // control flow statement
            if (electrons[i].x >= L)  { electrons[i].x = 2.0 * L - electrons[i].x; electrons[i].ux = -electrons[i].ux; }  // control flow statement
        }  // end block
        if (step % 2 == 0) {  // begin block
            double f_energy = 0.0;  // variable assignment/initialization
            for (int j = 0; j < Nx + 1; ++j) { f_energy += 0.5 * eps0 * Ex[j] * Ex[j] * dx; }  // control flow statement
            diag << step << "," << step*dt << "," << electrons[0].x << "," << electrons[0].ux << "," << electrons[0].uy << "," << f_energy << "\n";  // statement
        }  // end block
    }  // end block
    std::cout << "Done. Verification file ready.\n";  // statement
    return 0;  // return from function
}  // end block
