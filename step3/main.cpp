#include <iostream>  // include iostream
#include <vector>  // include vector
#include <cmath>  // include cmath
#include <fstream>  // include fstream
#include <omp.h>  // include omp.h
#include <algorithm>  // include algorithm
// (blank line)
struct Particle { double x, ux, uy, uz, w; };  // type definition
constexpr double PI = 3.14159265358979323846;  // variable assignment/initialization
constexpr double eps0 = 8.8541878128e-12, me = 9.1093837e-31, qe = -1.602176634e-19, c = 299792458.0;  // variable assignment/initialization
// (blank line)
int main() {  // begin block
    // Step 3 plan: Two-Stream Instability Growth Rate with OpenMP
    const int Nx = 512, ppc = 200, num_steps = 3000;  // variable assignment/initialization
    const double n_plasma = 1e19;  // variable assignment/initialization
    const double n_beam = 0.01 * n_plasma;  // variable assignment/initialization
    const double v_beam = 0.1 * c; // near-nonrelativistic for plan's analytical formula 
    const double gamma_b = 1.0 / std::sqrt(1.0 - v_beam * v_beam / (c * c));  // function declaration/call or variable initialization
// (blank line)
    const double omega_pe = std::sqrt(n_plasma * qe * qe / (eps0 * me));  // function declaration/call or variable initialization
    // Choose grid so that mode 1 is near resonance: k*v_beam ~ omega_pe
    const double k1 = omega_pe / v_beam;  // variable assignment/initialization
    const double L = 2.0 * PI / k1;  // variable assignment/initialization
    const double dx = L / Nx;  // variable assignment/initialization
    const double dt = 0.01 / omega_pe; // resolves plasma period 
// (blank line)
    std::vector<double> Ex(Nx, 0.0), rho(Nx, 0.0);  // function declaration/call or variable initialization
    std::vector<Particle> plasma_e(Nx * ppc);  // function declaration/call or variable initialization
    std::vector<Particle> beam_e(Nx * ppc);  // function declaration/call or variable initialization
// (blank line)
    const double w_p = (n_plasma * L) / plasma_e.size();  // function declaration/call or variable initialization
    const double w_b = (n_beam * L) / beam_e.size();  // function declaration/call or variable initialization
// (blank line)
    for (size_t i = 0; i < plasma_e.size(); ++i) {  // begin block
        plasma_e[i].x = ((i / ppc) + ((i % ppc) + 0.5) / ppc) * dx;  // variable assignment/initialization
        plasma_e[i].w = w_p; plasma_e[i].ux = 0.0; plasma_e[i].uy = 0.0; plasma_e[i].uz = 0.0;  // variable assignment/initialization
    }  // end block
    for (size_t i = 0; i < beam_e.size(); ++i) {  // begin block
        beam_e[i].x = ((i / ppc) + ((i % ppc) + 0.5) / ppc) * dx;  // variable assignment/initialization
        beam_e[i].w = w_b; beam_e[i].ux = gamma_b * v_beam; beam_e[i].uy = 0.0; beam_e[i].uz = 0.0;  // variable assignment/initialization
        // Small position perturbation to seed instability at mode 1
        beam_e[i].x += 1e-4 * dx * std::sin(2.0 * PI * beam_e[i].x / L);  // variable assignment/initialization
        if (beam_e[i].x >= L) beam_e[i].x -= L;  // control flow statement
        if (beam_e[i].x < 0.0) beam_e[i].x += L;  // control flow statement
    }  // end block
// (blank line)
    std::ofstream diag("two_stream_diagnostics.csv");  // function declaration/call or variable initialization
    diag << "step,time,E_k1_amp,instantaneous_gamma\n";  // statement
    double last_E_amp = 0.0;  // variable assignment/initialization
    double last_time = 0.0;  // variable assignment/initialization
// (blank line)
    for (int step = 1; step <= num_steps; ++step) {  // begin block
        // Charge deposition - serial for bitwise-deterministic results across OpenMP thread counts.
        std::fill(rho.begin(), rho.end(), 0.0);  // function declaration/call or variable initialization
        for (size_t i = 0; i < plasma_e.size(); ++i) {  // begin block
            double x_n = plasma_e[i].x / dx; int j = static_cast<int>(std::floor(x_n)); double xi = x_n - j;  // function declaration/call or variable initialization
            int j0 = (j % Nx + Nx) % Nx; int j1 = (j0 + 1) % Nx;  // function declaration/call or variable initialization
            rho[j0] += (qe * plasma_e[i].w * (1.0 - xi)) / dx;  // variable assignment/initialization
            rho[j1] += (qe * plasma_e[i].w * xi) / dx;  // variable assignment/initialization
        }  // end block
        for (size_t i = 0; i < beam_e.size(); ++i) {  // begin block
            double x_n = beam_e[i].x / dx; int j = static_cast<int>(std::floor(x_n)); double xi = x_n - j;  // function declaration/call or variable initialization
            int j0 = (j % Nx + Nx) % Nx; int j1 = (j0 + 1) % Nx;  // function declaration/call or variable initialization
            rho[j0] += (qe * beam_e[i].w * (1.0 - xi)) / dx;  // variable assignment/initialization
            rho[j1] += (qe * beam_e[i].w * xi) / dx;  // variable assignment/initialization
        }  // end block
// (blank line)
        // Neutralizing ion background (frozen)
        double total_e_charge = 0.0;  // variable assignment/initialization
        for (int j = 0; j < Nx; ++j) total_e_charge += rho[j];  // control flow statement
        double ion_bg = -total_e_charge / Nx;  // variable assignment/initialization
        for (int j = 0; j < Nx; ++j) rho[j] += ion_bg;  // control flow statement
// (blank line)
        // Spectral field solve (periodic)
        std::vector<double> rho_k_re(Nx, 0.0), rho_k_im(Nx, 0.0);  // function declaration/call or variable initialization
        for (int k = 0; k < Nx; ++k) {  // begin block
            double re = 0.0, im = 0.0;  // variable assignment/initialization
            for (int j = 0; j < Nx; ++j) {  // begin block
                double angle = 2.0 * PI * k * j / Nx;  // variable assignment/initialization
                re += rho[j] * std::cos(angle);  // variable assignment/initialization
                im -= rho[j] * std::sin(angle);  // variable assignment/initialization
            }  // end block
            rho_k_re[k] = re; rho_k_im[k] = im;  // variable assignment/initialization
        }  // end block
        std::vector<double> E_k_re(Nx, 0.0), E_k_im(Nx, 0.0);  // function declaration/call or variable initialization
        for (int k = 1; k < Nx; ++k) {  // begin block
            int k_signed = (k <= Nx/2) ? k : k - Nx;  // function declaration/call or variable initialization
            double kval = 2.0 * PI * k_signed / L;  // variable assignment/initialization
            // E_k = -i * rho_k / (eps0 * k)
            E_k_re[k] =  rho_k_im[k] / (eps0 * kval);  // variable assignment/initialization
            E_k_im[k] = -rho_k_re[k] / (eps0 * kval);  // variable assignment/initialization
        }  // end block
        if (Nx % 2 == 0) { E_k_re[Nx/2] = 0.0; E_k_im[Nx/2] = 0.0; }  // control flow statement
        std::fill(Ex.begin(), Ex.end(), 0.0);  // function declaration/call or variable initialization
        for (int j = 0; j < Nx; ++j) {  // begin block
            double val = 0.0;  // variable assignment/initialization
            for (int k = 0; k < Nx; ++k) {  // begin block
                double angle = 2.0 * PI * k * j / Nx;  // variable assignment/initialization
                val += E_k_re[k] * std::cos(angle) - E_k_im[k] * std::sin(angle);  // variable assignment/initialization
            }  // end block
            Ex[j] = val / Nx;  // variable assignment/initialization
        }  // end block
// (blank line)
        // Particle push - OpenMP parallelized (independent particle updates)
        #pragma omp parallel for schedule(static)  // preprocessor directive
        for (size_t i = 0; i < plasma_e.size(); ++i) {  // begin block
            double x_n = plasma_e[i].x / dx; int j = static_cast<int>(std::floor(x_n)); double xi = x_n - j;  // function declaration/call or variable initialization
            int j0 = (j % Nx + Nx) % Nx; int j1 = (j0 + 1) % Nx;  // function declaration/call or variable initialization
            double E = (1.0 - xi) * Ex[j0] + xi * Ex[j1];  // function declaration/call or variable initialization
            plasma_e[i].ux += (qe / me) * E * dt;  // variable assignment/initialization
            double g = std::sqrt(1.0 + (plasma_e[i].ux * plasma_e[i].ux) / (c * c));  // function declaration/call or variable initialization
            plasma_e[i].x  += (plasma_e[i].ux / g) * dt;  // variable assignment/initialization
            plasma_e[i].x = std::fmod(plasma_e[i].x, L);  // variable assignment/initialization
            if (plasma_e[i].x < 0.0) plasma_e[i].x += L;  // control flow statement
        }  // end block
// (blank line)
        #pragma omp parallel for schedule(static)  // preprocessor directive
        for (size_t i = 0; i < beam_e.size(); ++i) {  // begin block
            double x_n = beam_e[i].x / dx; int j = static_cast<int>(std::floor(x_n)); double xi = x_n - j;  // function declaration/call or variable initialization
            int j0 = (j % Nx + Nx) % Nx; int j1 = (j0 + 1) % Nx;  // function declaration/call or variable initialization
            double E = (1.0 - xi) * Ex[j0] + xi * Ex[j1];  // function declaration/call or variable initialization
            beam_e[i].ux += (qe / me) * E * dt;  // variable assignment/initialization
            double g = std::sqrt(1.0 + (beam_e[i].ux * beam_e[i].ux) / (c * c));  // function declaration/call or variable initialization
            beam_e[i].x  += (beam_e[i].ux / g) * dt;  // variable assignment/initialization
            beam_e[i].x = std::fmod(beam_e[i].x, L);  // variable assignment/initialization
            if (beam_e[i].x < 0.0) beam_e[i].x += L;  // control flow statement
        }  // end block
// (blank line)
        if (step % 5 == 0) {  // begin block
            double E_k1_re = 0.0, E_k1_im = 0.0;  // variable assignment/initialization
            for (int j = 0; j < Nx; ++j) {  // begin block
                double angle = 2.0 * PI * 1 * j / Nx;  // variable assignment/initialization
                E_k1_re += Ex[j] * std::cos(angle);  // variable assignment/initialization
                E_k1_im -= Ex[j] * std::sin(angle);  // variable assignment/initialization
            }  // end block
            double E_k1_amp = std::sqrt(E_k1_re*E_k1_re + E_k1_im*E_k1_im) / Nx;  // function declaration/call or variable initialization
            double current_time = step * dt;  // variable assignment/initialization
            double inst_gamma = 0.0;  // variable assignment/initialization
            if (last_E_amp > 0.0 && E_k1_amp > 0.0) {  // begin block
                inst_gamma = (std::log(E_k1_amp) - std::log(last_E_amp)) / (current_time - last_time);  // variable assignment/initialization
            }  // end block
            diag << step << "," << current_time << "," << E_k1_amp << "," << inst_gamma << "\n";  // statement
            last_E_amp = E_k1_amp;  // variable assignment/initialization
            last_time = current_time;  // variable assignment/initialization
        }  // end block
    }  // end block
    std::cout << "Step 3 complete. omega_pe=" << omega_pe << " rad/s\n";  // variable assignment/initialization
    return 0;  // return from function
}  // end block
