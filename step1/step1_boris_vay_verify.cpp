#include <iostream>  // include iostream
#include <vector>  // include vector
#include <cmath>  // include cmath
#include <fstream>  // include fstream
#include <iomanip>  // include iomanip
#include <algorithm>  // include algorithm
#include <numeric>  // include numeric
// (blank line)
constexpr double pi = 3.14159265358979323846;  // variable assignment/initialization
constexpr double eps0 = 8.8541878128e-12;  // variable assignment/initialization
constexpr double me   = 9.1093837015e-31;  // variable assignment/initialization
constexpr double qe   = -1.602176634e-19; // electron charge (negative) 
constexpr double c    = 299792458.0;  // variable assignment/initialization
// (blank line)
struct Particle {  // type definition
    double x, y, z;  // statement
    double ux, uy, uz; // gamma * v 
    double w;  // statement
};  // end block
// (blank line)
// Vay pusher (Vay, Phys. Plasmas 15, 056701, 2008)
// u = gamma * v, qdt2m = q * dt / (2 * m)
static inline void vay_push(Particle& p,  // code line
                            double Ex, double Ey, double Ez,  // code line
                            double Bx, double By, double Bz,  // code line
                            double qdt2m) {  // begin block
    // Half electric acceleration
    double ux_minus = p.ux + qdt2m * Ex;  // variable assignment/initialization
    double uy_minus = p.uy + qdt2m * Ey;  // variable assignment/initialization
    double uz_minus = p.uz + qdt2m * Ez;  // variable assignment/initialization
// (blank line)
    // Magnetic rotation (Vay)
    double gamma_minus_sq =  // code line
        1.0 + (ux_minus * ux_minus + uy_minus * uy_minus + uz_minus * uz_minus) / (c * c);  // statement
// (blank line)
    double tx = qdt2m * Bx;  // variable assignment/initialization
    double ty = qdt2m * By;  // variable assignment/initialization
    double tz = qdt2m * Bz;  // variable assignment/initialization
    double tau_sq = tx * tx + ty * ty + tz * tz;  // variable assignment/initialization
// (blank line)
    double u_star = (ux_minus * tx + uy_minus * ty + uz_minus * tz) / c; // u^- . tau / c 
    double sigma = gamma_minus_sq - tau_sq;  // variable assignment/initialization
    double gamma_plus = std::sqrt(  // code line
        0.5 * (sigma + std::sqrt(sigma * sigma + 4.0 * (tau_sq + u_star * u_star))));  // statement
// (blank line)
    double inv_gamma_plus = 1.0 / gamma_plus;  // variable assignment/initialization
    double sx = tx * inv_gamma_plus;  // variable assignment/initialization
    double sy = ty * inv_gamma_plus;  // variable assignment/initialization
    double sz = tz * inv_gamma_plus;  // variable assignment/initialization
    double s_sq = sx * sx + sy * sy + sz * sz;  // variable assignment/initialization
    double coeff = 2.0 / (1.0 + s_sq);  // function declaration/call or variable initialization
// (blank line)
    // u' = u^- + u^- x s
    double up_x = ux_minus + (uy_minus * sz - uz_minus * sy);  // function declaration/call or variable initialization
    double up_y = uy_minus + (uz_minus * sx - ux_minus * sz);  // function declaration/call or variable initialization
    double up_z = uz_minus + (ux_minus * sy - uy_minus * sx);  // function declaration/call or variable initialization
// (blank line)
    // u^+ = u^- + coeff * (u' x s)
    double uplus_x = ux_minus + coeff * (up_y * sz - up_z * sy);  // function declaration/call or variable initialization
    double uplus_y = uy_minus + coeff * (up_z * sx - up_x * sz);  // function declaration/call or variable initialization
    double uplus_z = uz_minus + coeff * (up_x * sy - up_y * sx);  // function declaration/call or variable initialization
// (blank line)
    // Half electric acceleration
    p.ux = uplus_x + qdt2m * Ex;  // variable assignment/initialization
    p.uy = uplus_y + qdt2m * Ey;  // variable assignment/initialization
    p.uz = uplus_z + qdt2m * Ez;  // variable assignment/initialization
}  // end block
// (blank line)
static inline double gamma_of(const Particle& p) {  // begin block
    double u2 = p.ux * p.ux + p.uy * p.uy + p.uz * p.uz;  // variable assignment/initialization
    return std::sqrt(1.0 + u2 / (c * c));  // return from function
}  // end block
// (blank line)
int main() {  // begin block
    // Step 1 physics check from plan:
    // single electron, gamma = 2.0, perpendicular velocity, uniform B0 = B0*z, E=0.
    const double gamma0 = 2.0;  // variable assignment/initialization
    const double B0z = 1.0; // Tesla, uniform static magnetic field 
    const double Ex = 0.0, Ey = 0.0, Ez = 0.0;  // variable assignment/initialization
    const double Bx = 0.0, By = 0.0, Bz = B0z;  // variable assignment/initialization
// (blank line)
    const double v_perp = c * std::sqrt(1.0 - 1.0 / (gamma0 * gamma0)); // |v| for gamma=2 
    const double u_perp = gamma0 * v_perp; // gamma*v 
// (blank line)
    Particle p{};  // statement
    p.x = 0.0;  // variable assignment/initialization
    p.y = 0.0;  // variable assignment/initialization
    p.z = 0.0;  // variable assignment/initialization
    p.ux = u_perp; // perpendicular to Bz 
    p.uy = 0.0;  // variable assignment/initialization
    p.uz = 0.0;  // variable assignment/initialization
    p.w = 1.0;  // variable assignment/initialization
// (blank line)
    // Analytical values
    const double RL_analytic = (gamma0 * me * v_perp) / (std::abs(qe) * B0z);  // function declaration/call or variable initialization
    const double omega_c_analytic = (std::abs(qe) * B0z) / (gamma0 * me);  // function declaration/call or variable initialization
    const double Tc_analytic = 2.0 * pi / omega_c_analytic;  // variable assignment/initialization
// (blank line)
    // Time integration: 10^4 steps, resolve cyclotron motion well.
    const int num_steps = 10000;  // variable assignment/initialization
    const double dt = Tc_analytic / 200.0; // 200 steps per cyclotron period 
    const double qdt2m = qe * dt / (2.0 * me);  // function declaration/call or variable initialization
// (blank line)
    std::ofstream traj("step1_trajectory.csv");  // function declaration/call or variable initialization
    traj << std::setprecision(17);  // statement
    traj << "step,time,x,y,z,ux,uy,uz,gamma\n";  // statement
// (blank line)
    std::vector<double> xs;  // statement
    std::vector<double> ys;  // statement
    xs.reserve(num_steps + 1);  // statement
    ys.reserve(num_steps + 1);  // statement
// (blank line)
    double gamma_min = 1e300;  // variable assignment/initialization
    double gamma_max = -1e300;  // variable assignment/initialization
    double gamma_sum = 0.0;  // variable assignment/initialization
// (blank line)
    // Record initial state
    auto record = [&](int step) {  // begin block
        double g = gamma_of(p);  // function declaration/call or variable initialization
        gamma_min = std::min(gamma_min, g);  // variable assignment/initialization
        gamma_max = std::max(gamma_max, g);  // variable assignment/initialization
        gamma_sum += g;  // variable assignment/initialization
        xs.push_back(p.x);  // statement
        ys.push_back(p.y);  // statement
        traj << step << ',' << step * dt << ',' << p.x << ',' << p.y << ',' << p.z << ','  // code line
             << p.ux << ',' << p.uy << ',' << p.uz << ',' << g << '\n';  // statement
    };  // end block
// (blank line)
    record(0);  // statement
    for (int step = 1; step <= num_steps; ++step) {  // begin block
        vay_push(p, Ex, Ey, Ez, Bx, By, Bz, qdt2m);  // statement
        double g = gamma_of(p);  // function declaration/call or variable initialization
        // position advance: dx = v * dt = (u/gamma) * dt
        p.x += (p.ux / g) * dt;  // variable assignment/initialization
        p.y += (p.uy / g) * dt;  // variable assignment/initialization
        p.z += (p.uz / g) * dt;  // variable assignment/initialization
        record(step);  // statement
    }  // end block
    traj.close();  // statement
// (blank line)
    // Radius from trajectory extent (circle should be closed over many periods)
    double xmin = *std::min_element(xs.begin(), xs.end());  // function declaration/call or variable initialization
    double xmax = *std::max_element(xs.begin(), xs.end());  // function declaration/call or variable initialization
    double ymin = *std::min_element(ys.begin(), ys.end());  // function declaration/call or variable initialization
    double ymax = *std::max_element(ys.begin(), ys.end());  // function declaration/call or variable initialization
    double cx_est = 0.5 * (xmin + xmax);  // function declaration/call or variable initialization
    double cy_est = 0.5 * (ymin + ymax);  // function declaration/call or variable initialization
    double Rx = 0.5 * (xmax - xmin);  // function declaration/call or variable initialization
    double Ry = 0.5 * (ymax - ymin);  // function declaration/call or variable initialization
    double RL_measured_extent = 0.5 * (Rx + Ry);  // function declaration/call or variable initialization
// (blank line)
    // Radius from momentum magnitude (should be invariant in uniform B, E=0)
    double u_perp_measured = std::sqrt(p.ux * p.ux + p.uy * p.uy);  // function declaration/call or variable initialization
    double RL_momentum = (me * u_perp_measured) / (std::abs(qe) * B0z);  // function declaration/call or variable initialization
// (blank line)
    // Omega_c from phase unwrapping around estimated center
    std::vector<double> phases;  // statement
    phases.reserve(xs.size());  // statement
    for (size_t i = 0; i < xs.size(); ++i) {  // begin block
        phases.push_back(std::atan2(ys[i] - cy_est, xs[i] - cx_est));  // statement
    }  // end block
    // unwrap
    double unwrapped = phases[0];  // variable assignment/initialization
    double prev = phases[0];  // variable assignment/initialization
    double phase_start = unwrapped;  // variable assignment/initialization
    for (size_t i = 1; i < phases.size(); ++i) {  // begin block
        double d = phases[i] - prev;  // variable assignment/initialization
        if (d > pi) d -= 2.0 * pi;  // control flow statement
        if (d < -pi) d += 2.0 * pi;  // control flow statement
        unwrapped += d;  // variable assignment/initialization
        prev = phases[i];  // variable assignment/initialization
    }  // end block
    double phase_end = unwrapped;  // variable assignment/initialization
    double total_time = num_steps * dt;  // variable assignment/initialization
    double omega_c_measured = std::abs(phase_end - phase_start) / total_time;  // function declaration/call or variable initialization
// (blank line)
    // Energy/gamma invariance
    double gamma_final = gamma_of(p);  // function declaration/call or variable initialization
    double gamma_drift = std::abs(gamma_final - gamma0) / gamma0;  // function declaration/call or variable initialization
    double gamma_range = (gamma_max - gamma_min) / gamma0;  // function declaration/call or variable initialization
// (blank line)
    // Errors
    double RL_err_extent = std::abs(RL_measured_extent - RL_analytic) / RL_analytic * 100.0;  // function declaration/call or variable initialization
    double RL_err_momentum = std::abs(RL_momentum - RL_analytic) / RL_analytic * 100.0;  // function declaration/call or variable initialization
    double omega_err = std::abs(omega_c_measured - omega_c_analytic) / omega_c_analytic * 100.0;  // function declaration/call or variable initialization
// (blank line)
    bool pass_RL = (RL_err_extent < 0.1);  // function declaration/call or variable initialization
    bool pass_omega = (omega_err < 0.1);  // function declaration/call or variable initialization
    bool pass_gamma = (gamma_range <= 1e-12) && (gamma_drift <= 1e-12);  // function declaration/call or variable initialization
// (blank line)
    std::cout << std::setprecision(12);  // function declaration/call or variable initialization
    std::cout << "============================================================\n";  // variable assignment/initialization
    std::cout << " Step 1: Relativistic Boris-Vay Pusher Verification\n";  // statement
    std::cout << "============================================================\n";  // variable assignment/initialization
    std::cout << "Parameters:\n";  // statement
    std::cout << "  gamma0                = " << gamma0 << "\n";  // variable assignment/initialization
    std::cout << "  B0z (T)               = " << B0z << "\n";  // function declaration/call or variable initialization
    std::cout << "  dt (s)                = " << dt << "\n";  // function declaration/call or variable initialization
    std::cout << "  steps                 = " << num_steps << "\n";  // variable assignment/initialization
    std::cout << "  total time / Tc       = " << total_time / Tc_analytic << "\n";  // variable assignment/initialization
    std::cout << "\nAnalytical:\n";  // statement
    std::cout << "  R_L (m)               = " << RL_analytic << "\n";  // function declaration/call or variable initialization
    std::cout << "  omega_c (rad/s)       = " << omega_c_analytic << "\n";  // function declaration/call or variable initialization
    std::cout << "  Tc (s)                = " << Tc_analytic << "\n";  // function declaration/call or variable initialization
    std::cout << "\nMeasured:\n";  // statement
    std::cout << "  R_L extent (m)        = " << RL_measured_extent << "  (err=" << RL_err_extent << "%)\n";  // function declaration/call or variable initialization
    std::cout << "  R_L momentum (m)      = " << RL_momentum << "  (err=" << RL_err_momentum << "%)\n";  // function declaration/call or variable initialization
    std::cout << "  omega_c (rad/s)       = " << omega_c_measured << "  (err=" << omega_err << "%)\n";  // function declaration/call or variable initialization
    std::cout << "  gamma initial         = " << gamma0 << "\n";  // variable assignment/initialization
    std::cout << "  gamma final           = " << gamma_final << "\n";  // variable assignment/initialization
    std::cout << "  gamma min             = " << gamma_min << "\n";  // variable assignment/initialization
    std::cout << "  gamma max             = " << gamma_max << "\n";  // variable assignment/initialization
    std::cout << "  gamma drift (rel)     = " << gamma_drift << "\n";  // function declaration/call or variable initialization
    std::cout << "  gamma range (rel)     = " << gamma_range << "\n";  // function declaration/call or variable initialization
    std::cout << "\nVerification:\n";  // statement
    std::cout << "  R_L within 0.1%       : " << (pass_RL ? "PASS" : "FAIL") << "\n";  // function declaration/call or variable initialization
    std::cout << "  omega_c within 0.1%   : " << (pass_omega ? "PASS" : "FAIL") << "\n";  // function declaration/call or variable initialization
    std::cout << "  gamma invariant       : " << (pass_gamma ? "PASS" : "FAIL") << "\n";  // function declaration/call or variable initialization
    std::cout << "  OVERALL               : " << ((pass_RL && pass_omega && pass_gamma) ? "PASS" : "FAIL") << "\n";  // function declaration/call or variable initialization
    std::cout << "============================================================\n";  // variable assignment/initialization
    std::cout << "Trajectory written to step1_trajectory.csv\n";  // statement
// (blank line)
    return (pass_RL && pass_omega && pass_gamma) ? 0 : 1;  // return from function
}  // end block
