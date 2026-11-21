#include <iostream>  // include iostream
#include <vector>  // include vector
#include <cmath>  // include cmath
#include <fstream>  // include fstream
#include <iomanip>  // include iomanip
#include <algorithm>  // include algorithm
#include <random>  // include random
#include <numeric>  // include numeric
// (blank line)
constexpr double pi = 3.14159265358979323846;  // variable assignment/initialization
constexpr double eps0 = 8.8541878128e-12;  // variable assignment/initialization
constexpr double me   = 9.1093837015e-31;  // variable assignment/initialization
constexpr double qe   = -1.602176634e-19; // electron charge (negative) 
constexpr double c    = 299792458.0;  // variable assignment/initialization
// (blank line)
struct Particle {  // type definition
    double x;  // statement
    double ux, uy, uz; // gamma*v 
    double w;  // statement
};  // end block
// (blank line)
static inline double gamma_of(const Particle& p) {  // begin block
    double u2 = p.ux * p.ux + p.uy * p.uy + p.uz * p.uz;  // variable assignment/initialization
    return std::sqrt(1.0 + u2 / (c * c));  // return from function
}  // end block
// (blank line)
// Cooley-Tukey radix-2 FFT (in-place)
static void fft(std::vector<double>& re, std::vector<double>& im, bool inverse) {  // begin block
    const int n = static_cast<int>(re.size());  // function declaration/call or variable initialization
    for (int i = 1, j = 0; i < n; ++i) {  // begin block
        int bit = n >> 1;  // variable assignment/initialization
        for (; j & bit; bit >>= 1) j ^= bit;  // control flow statement
        j ^= bit;  // variable assignment/initialization
        if (i < j) { std::swap(re[i], re[j]); std::swap(im[i], im[j]); }  // control flow statement
    }  // end block
    for (int len = 2; len <= n; len <<= 1) {  // begin block
        double ang = 2.0 * pi / len * (inverse ? 1.0 : -1.0);  // function declaration/call or variable initialization
        double wlen_r = std::cos(ang), wlen_i = std::sin(ang);  // function declaration/call or variable initialization
        for (int i = 0; i < n; i += len) {  // begin block
            double wr = 1.0, wi = 0.0;  // variable assignment/initialization
            for (int j = 0; j < len / 2; ++j) {  // begin block
                int u = i + j;  // variable assignment/initialization
                int v = i + j + len / 2;  // variable assignment/initialization
                double vr = re[v] * wr - im[v] * wi;  // variable assignment/initialization
                double vi = re[v] * wi + im[v] * wr;  // variable assignment/initialization
                re[v] = re[u] - vr; im[v] = im[u] - vi;  // variable assignment/initialization
                re[u] += vr;       im[u] += vi;  // variable assignment/initialization
                double nwr = wr * wlen_r - wi * wlen_i;  // variable assignment/initialization
                wi = wr * wlen_i + wi * wlen_r;  // variable assignment/initialization
                wr = nwr;  // variable assignment/initialization
            }  // end block
        }  // end block
    }  // end block
    if (inverse) {  // begin block
        for (int i = 0; i < n; ++i) { re[i] /= n; im[i] /= n; }  // control flow statement
    }  // end block
}  // end block
// (blank line)
// Quadratic (TSC) shape function weights for distance xi = x/dx - j, |xi| < 1.5
static inline void tsc_weights(double xi, double& wm, double& w0, double& wp) {  // begin block
    // Cloud-in-cell/TSC: S(xi)=3/4-xi^2 for |xi|<1/2, 1/2*(3/2-|xi|)^2 for 1/2<=|xi|<3/2
    // with xi relative to nearest grid point j.
    wm = 0.5 * (0.5 - xi) * (0.5 - xi);  // variable assignment/initialization
    w0 = 0.75 - xi * xi;  // variable assignment/initialization
    wp = 0.5 * (0.5 + xi) * (0.5 + xi);  // variable assignment/initialization
}  // end block
// (blank line)
int main() {  // begin block
    // --- Parameters chosen for clear omega_pe signal and manageable runtime ---
    const int Nx = 256;                 // power of two for FFT 
    const int ppc = 200;                // particles per cell 
    const double n0 = 1e16;             // m^-3, cold plasma density 
    const double Te_eV = 1.0;           // small thermal temperature (eV) 
    const double delta = 0.01;          // perturbation amplitude 
    const int perturb_mode = 1;         // k = 2*pi*mode/L 
// (blank line)
    const double omega_pe = std::sqrt(n0 * qe * qe / (eps0 * me));  // function declaration/call or variable initialization
    const double Tpe = 2.0 * pi / omega_pe;  // variable assignment/initialization
// (blank line)
    // Resolve Debye length and plasma period; domain several Debye lengths.
    const double vth = std::sqrt(2.0 * Te_eV * std::abs(qe) / me); // thermal speed 
    const double lambda_D = std::sqrt(eps0 * Te_eV * std::abs(qe) / (n0 * qe * qe));  // function declaration/call or variable initialization
    double dx = std::max(lambda_D, 1e-4); // ensure reasonable cell size 
    double L = Nx * dx;  // variable assignment/initialization
    double k_perturb = 2.0 * pi * perturb_mode / L;  // variable assignment/initialization
// (blank line)
    // dt resolves omega_pe well and satisfies CFL-ish constraint for thermal motion
    double dt = 0.02 / omega_pe;  // variable assignment/initialization
    // Ensure dt small enough for particle not to cross many cells thermally
    double dt_cfl = 0.2 * dx / std::max(vth, 1e3);  // function declaration/call or variable initialization
    if (dt > dt_cfl) dt = dt_cfl;  // control flow statement
    // Make total time 1000 omega_pe^-1 as requested
    const double total_time = 1000.0 / omega_pe;  // variable assignment/initialization
    const int num_steps = static_cast<int>(std::ceil(total_time / dt));  // function declaration/call or variable initialization
    // Redefine dt to hit total_time exactly in integer steps? keep simple
    dt = total_time / num_steps;  // variable assignment/initialization
// (blank line)
    std::cout << std::setprecision(10);  // function declaration/call or variable initialization
    std::cout << "Step 2: Plasma oscillation / high-order deposition verification\n";  // statement
    std::cout << "  Nx=" << Nx << " ppc=" << ppc << " n0=" << n0 << " Te_eV=" << Te_eV << "\n";  // variable assignment/initialization
    std::cout << "  dx=" << dx << " L=" << L << " dt=" << dt << " steps=" << num_steps << "\n";  // variable assignment/initialization
    std::cout << "  omega_pe=" << omega_pe << " rad/s, Tpe=" << Tpe << " s\n";  // variable assignment/initialization
    std::cout << "  lambda_D=" << lambda_D << " m, vth=" << vth << " m/s\n";  // variable assignment/initialization
    std::cout << "  perturb_mode=" << perturb_mode << " delta=" << delta << "\n";  // variable assignment/initialization
// (blank line)
    // --- Particle initialization ---
    const int Np = Nx * ppc;  // variable assignment/initialization
    std::vector<Particle> electrons(Np);  // function declaration/call or variable initialization
    const double macro_weight = n0 * L / Np; // represents number of real electrons 
    std::mt19937 rng(12345);  // function declaration/call or variable initialization
    std::normal_distribution<double> normal(0.0, 1.0);  // function declaration/call or variable initialization
// (blank line)
    for (int i = 0; i < Np; ++i) {  // begin block
        double x0 = (i + 0.5) * L / Np;  // function declaration/call or variable initialization
        // Position perturbation to create density perturbation:
        // x = x0 + (delta/k) * sin(k*x0) => n_e ~ n0*(1 - delta*cos(kx)) approx; use exact target form by shifting
        double x = x0 + (delta / k_perturb) * std::sin(k_perturb * x0);  // function declaration/call or variable initialization
        // wrap
        x = std::fmod(x, L);  // variable assignment/initialization
        if (x < 0) x += L;  // control flow statement
        electrons[i].x = x;  // variable assignment/initialization
        electrons[i].w = macro_weight;  // variable assignment/initialization
        // Maxwellian in ux,uy,uz (non-relativistic for 1 eV, but store u=gamma*v)
        double vx = normal(rng) * std::sqrt(Te_eV * std::abs(qe) / me);  // function declaration/call or variable initialization
        double vy = normal(rng) * std::sqrt(Te_eV * std::abs(qe) / me);  // function declaration/call or variable initialization
        double vz = normal(rng) * std::sqrt(Te_eV * std::abs(qe) / me);  // function declaration/call or variable initialization
        double v2 = vx*vx + vy*vy + vz*vz;  // variable assignment/initialization
        double g = 1.0 / std::sqrt(std::max(1e-30, 1.0 - v2/(c*c)));  // function declaration/call or variable initialization
        electrons[i].ux = g * vx;  // variable assignment/initialization
        electrons[i].uy = g * vy;  // variable assignment/initialization
        electrons[i].uz = g * vz;  // variable assignment/initialization
    }  // end block
// (blank line)
    // --- Grid arrays ---
    std::vector<double> rho(Nx, 0.0);  // function declaration/call or variable initialization
    std::vector<double> Ex(Nx, 0.0);  // function declaration/call or variable initialization
// (blank line)
    // Frozen ion background: charge density +n0*|qe| = -n0*qe
    const double rho_ion = -n0 * qe;  // variable assignment/initialization
// (blank line)
    // --- Diagnostics ---
    std::ofstream diag("step2_diagnostics.csv");  // function declaration/call or variable initialization
    diag << std::setprecision(17);  // statement
    diag << "step,time,field_energy,E_mode_re,E_mode_im,kinetic_energy,total_energy\n";  // statement
// (blank line)
    std::vector<double> sample_t;  // statement
    std::vector<double> sample_Ere;  // statement
    std::vector<double> sample_Eim;  // statement
    std::vector<double> sample_W;  // statement
    std::vector<double> sample_KE;  // statement
// (blank line)
    // Sampling: enough points for FFT, target ~10 points per plasma period
    int sample_every = std::max(1, static_cast<int>((Tpe / 20.0) / dt));  // function declaration/call or variable initialization
    std::cout << "  sample_every=" << sample_every << " steps\n";  // variable assignment/initialization
// (blank line)
    auto deposit = [&](void) {  // begin block
        std::fill(rho.begin(), rho.end(), rho_ion);  // function declaration/call or variable initialization
        for (const auto& p : electrons) {  // begin block
            double x_norm = p.x / dx; // in [0,Nx) 
            int j = static_cast<int>(std::floor(x_norm + 0.5)); // nearest grid point 
            double xi = x_norm - j;  // variable assignment/initialization
            // TSC weights around j: j-1,j,j+1
            double wm, w0, wp;  // statement
            tsc_weights(xi, wm, w0, wp);  // statement
            int jm = (j - 1 + Nx) % Nx;  // function declaration/call or variable initialization
            int j0 = (j + Nx) % Nx;  // function declaration/call or variable initialization
            int jp = (j + 1) % Nx;  // function declaration/call or variable initialization
            double qdens = qe * p.w / dx;  // variable assignment/initialization
            rho[jm] += qdens * wm;  // variable assignment/initialization
            rho[j0] += qdens * w0;  // variable assignment/initialization
            rho[jp] += qdens * wp;  // variable assignment/initialization
        }  // end block
    };  // end block
// (blank line)
    auto solve_field = [&](void) {  // begin block
        // Spectral Poisson solve: -i*k*phi_k = rho_k/eps0? Better: dE/dx=rho/eps0, E_k = -i*rho_k/(eps0*k)
        std::vector<double> re(Nx, 0.0), im(Nx, 0.0);  // function declaration/call or variable initialization
        for (int i = 0; i < Nx; ++i) re[i] = rho[i];  // control flow statement
        fft(re, im, false);  // statement
        // k=0: set to 0 (neutral)
        re[0] = 0.0; im[0] = 0.0;  // variable assignment/initialization
        for (int k = 1; k < Nx; ++k) {  // begin block
            int kidx = (k <= Nx/2) ? k : k - Nx;  // function declaration/call or variable initialization
            double kval = 2.0 * pi * kidx / L;  // variable assignment/initialization
            // E_k = -i * rho_k / (eps0 * kval)
            double rho_r = re[k], rho_i = im[k];  // variable assignment/initialization
            re[k] =  rho_i / (eps0 * kval);  // variable assignment/initialization
            im[k] = -rho_r / (eps0 * kval);  // variable assignment/initialization
        }  // end block
        // For even Nx, Nyquist mode - set to 0
        if (Nx % 2 == 0) { re[Nx/2] = 0.0; im[Nx/2] = 0.0; }  // control flow statement
        fft(re, im, true);  // statement
        for (int i = 0; i < Nx; ++i) Ex[i] = re[i];  // control flow statement
    };  // end block
// (blank line)
    auto gather_Ex = [&](double x) -> double {  // begin block
        double x_norm = x / dx;  // variable assignment/initialization
        int j = static_cast<int>(std::floor(x_norm + 0.5));  // function declaration/call or variable initialization
        double xi = x_norm - j;  // variable assignment/initialization
        double wm, w0, wp;  // statement
        tsc_weights(xi, wm, w0, wp);  // statement
        int jm = (j - 1 + Nx) % Nx;  // function declaration/call or variable initialization
        int j0 = (j + Nx) % Nx;  // function declaration/call or variable initialization
        int jp = (j + 1) % Nx;  // function declaration/call or variable initialization
        return wm * Ex[jm] + w0 * Ex[j0] + wp * Ex[jp];  // return from function
    };  // end block
// (blank line)
    // Initial field
    deposit();  // statement
    solve_field();  // statement
// (blank line)
    for (int step = 0; step <= num_steps; ++step) {  // begin block
        // Diagnostics
        if (step % sample_every == 0 || step == num_steps) {  // begin block
            double W = 0.0;  // variable assignment/initialization
            for (int j = 0; j < Nx; ++j) W += 0.5 * eps0 * Ex[j] * Ex[j] * dx;  // control flow statement
            double KE = 0.0;  // variable assignment/initialization
            for (const auto& p : electrons) {  // begin block
                double g = gamma_of(p);  // function declaration/call or variable initialization
                KE += (g - 1.0) * me * c * c * p.w;  // variable assignment/initialization
            }  // end block
            // Dominant mode amplitude (mode 1) via DFT of Ex
            double Er = 0.0, Ei = 0.0;  // variable assignment/initialization
            for (int j = 0; j < Nx; ++j) {  // begin block
                double angle = -2.0 * pi * perturb_mode * j / Nx;  // variable assignment/initialization
                Er += Ex[j] * std::cos(angle);  // variable assignment/initialization
                Ei += Ex[j] * std::sin(angle);  // variable assignment/initialization
            }  // end block
            Er /= Nx; Ei /= Nx;  // variable assignment/initialization
            double t = step * dt;  // variable assignment/initialization
            diag << step << ',' << t << ',' << W << ',' << Er << ',' << Ei << ',' << KE << ',' << (W + KE) << '\n';  // statement
            sample_t.push_back(t);  // statement
            sample_W.push_back(W);  // statement
            sample_Ere.push_back(Er);  // statement
            sample_Eim.push_back(Ei);  // statement
            sample_KE.push_back(KE);  // statement
        }  // end block
// (blank line)
        if (step == num_steps) break;  // control flow statement
// (blank line)
        // Boris push (relativistic, E only, B=0) - leapfrog style: we have x^n, u^{n-1/2}? 
        // For simplicity use Boris/Vay-like with E-only: u^{n+1/2}=u^{n-1/2}+q/m*E^n*dt, x^{n+1}=x^n+u/g*dt
        // Initialize by half-step back at step 0
        if (step == 0) {  // begin block
            // Backward half-kick so leapfrog is centered
            for (auto& p : electrons) {  // begin block
                double E = gather_Ex(p.x);  // function declaration/call or variable initialization
                p.ux -= 0.5 * (qe / me) * E * dt;  // variable assignment/initialization
            }  // end block
        }  // end block
// (blank line)
        // Full kick
        for (auto& p : electrons) {  // begin block
            double E = gather_Ex(p.x);  // function declaration/call or variable initialization
            p.ux += (qe / me) * E * dt; // uy,uz unchanged (E_x only) 
        }  // end block
        // Drift
        for (auto& p : electrons) {  // begin block
            double g = gamma_of(p);  // function declaration/call or variable initialization
            p.x += (p.ux / g) * dt;  // variable assignment/initialization
            // periodic
            p.x = std::fmod(p.x, L);  // variable assignment/initialization
            if (p.x < 0) p.x += L;  // control flow statement
        }  // end block
// (blank line)
        deposit();  // statement
        solve_field();  // statement
    }  // end block
    diag.close();  // statement
// (blank line)
    // --- Post-processing: FFT of E_mode real time history ---
    // Use sample_Ere as E-field time history for dominant mode
    int Nfft = 1;  // variable assignment/initialization
    while (Nfft < static_cast<int>(sample_Ere.size())) Nfft <<= 1;  // control flow statement
    std::vector<double> fr(Nfft, 0.0), fi(Nfft, 0.0);  // function declaration/call or variable initialization
    for (size_t i = 0; i < sample_Ere.size(); ++i) fr[i] = sample_Ere[i];  // control flow statement
    // remove DC
    double mean = std::accumulate(sample_Ere.begin(), sample_Ere.end(), 0.0) / sample_Ere.size();  // function declaration/call or variable initialization
    for (size_t i = 0; i < sample_Ere.size(); ++i) fr[i] = sample_Ere[i] - mean;  // control flow statement
    fft(fr, fi, false);  // statement
// (blank line)
    double sample_dt = (sample_t.size() > 1) ? (sample_t.back() - sample_t.front()) / (sample_t.size() - 1) : dt;  // function declaration/call or variable initialization
    double df = 1.0 / (Nfft * sample_dt);  // function declaration/call or variable initialization
    int peak_bin = 1;  // variable assignment/initialization
    double peak_amp = 0.0;  // variable assignment/initialization
    for (int k = 1; k < Nfft / 2; ++k) {  // begin block
        double amp = std::sqrt(fr[k]*fr[k] + fi[k]*fi[k]);  // function declaration/call or variable initialization
        if (amp > peak_amp) { peak_amp = amp; peak_bin = k; }  // control flow statement
    }  // end block
    double f_peak = peak_bin * df;  // variable assignment/initialization
    double omega_peak = 2.0 * pi * f_peak;  // variable assignment/initialization
    double omega_err_pct = std::abs(omega_peak - omega_pe) / omega_pe * 100.0;  // function declaration/call or variable initialization
// (blank line)
    // Also check energy oscillation frequency (should be 2*omega_pe)
    int NfftW = 1;  // variable assignment/initialization
    while (NfftW < static_cast<int>(sample_W.size())) NfftW <<= 1;  // control flow statement
    std::vector<double> Wr(NfftW, 0.0), Wi(NfftW, 0.0);  // function declaration/call or variable initialization
    double meanW = std::accumulate(sample_W.begin(), sample_W.end(), 0.0) / sample_W.size();  // function declaration/call or variable initialization
    for (size_t i = 0; i < sample_W.size(); ++i) Wr[i] = sample_W[i] - meanW;  // control flow statement
    fft(Wr, Wi, false);  // statement
    int peak_bin_W = 1; double peak_amp_W = 0.0;  // variable assignment/initialization
    for (int k = 1; k < NfftW / 2; ++k) {  // begin block
        double amp = std::sqrt(Wr[k]*Wr[k] + Wi[k]*Wi[k]);  // function declaration/call or variable initialization
        if (amp > peak_amp_W) { peak_amp_W = amp; peak_bin_W = k; }  // control flow statement
    }  // end block
    double omega_peak_W = 2.0 * pi * peak_bin_W * df;  // variable assignment/initialization
// (blank line)
    // Heating check: compare KE over early and late windows
    size_t nS = sample_KE.size();  // function declaration/call or variable initialization
    size_t w = std::max<size_t>(1, nS / 10);  // function declaration/call or variable initialization
    double KE_early = std::accumulate(sample_KE.begin(), sample_KE.begin() + w, 0.0) / w;  // function declaration/call or variable initialization
    double KE_late  = std::accumulate(sample_KE.end() - w, sample_KE.end(), 0.0) / w;  // function declaration/call or variable initialization
    double KE_change_pct = std::abs(KE_late - KE_early) / std::max(1e-30, KE_early) * 100.0;  // function declaration/call or variable initialization
    // Total energy drift
    double TE0 = sample_W.front() + sample_KE.front();  // function declaration/call or variable initialization
    double TE1 = sample_W.back() + sample_KE.back();  // function declaration/call or variable initialization
    double TE_drift_pct = std::abs(TE1 - TE0) / std::max(1e-30, std::abs(TE0)) * 100.0;  // function declaration/call or variable initialization
// (blank line)
    bool pass_omega = omega_err_pct < 5.0; // allow some spectral resolution tolerance 
    bool pass_heating = KE_change_pct < 10.0; // flat over 1000 omega_pe^-1 for this setup 
// (blank line)
    std::cout << "\n--- Verification ---\n";  // statement
    std::cout << "  omega_pe (analytical) = " << omega_pe << " rad/s\n";  // function declaration/call or variable initialization
    std::cout << "  FFT peak (E_mode)     = " << omega_peak << " rad/s  (err=" << omega_err_pct << "%)\n";  // function declaration/call or variable initialization
    std::cout << "  FFT peak (W_E)        = " << omega_peak_W << " rad/s  (~2*omega_pe=" << 2.0*omega_pe << ")\n";  // function declaration/call or variable initialization
    std::cout << "  KE early              = " << KE_early << " J\n";  // variable assignment/initialization
    std::cout << "  KE late               = " << KE_late << " J\n";  // variable assignment/initialization
    std::cout << "  KE change             = " << KE_change_pct << "%\n";  // variable assignment/initialization
    std::cout << "  Total energy drift    = " << TE_drift_pct << "%\n";  // variable assignment/initialization
    std::cout << "  omega_pe peak <5% err : " << (pass_omega ? "PASS" : "FAIL") << "\n";  // function declaration/call or variable initialization
    std::cout << "  heating flat <10%     : " << (pass_heating ? "PASS" : "FAIL") << "\n";  // function declaration/call or variable initialization
    std::cout << "  OVERALL               : " << ((pass_omega && pass_heating) ? "PASS" : "FAIL") << "\n";  // function declaration/call or variable initialization
// (blank line)
    std::ofstream summary("step2_verification.txt");  // function declaration/call or variable initialization
    summary << std::setprecision(12);  // statement
    summary << "omega_pe_analytic=" << omega_pe << "\n";  // variable assignment/initialization
    summary << "omega_peak_Emode=" << omega_peak << "\n";  // variable assignment/initialization
    summary << "omega_err_pct=" << omega_err_pct << "\n";  // variable assignment/initialization
    summary << "omega_peak_WE=" << omega_peak_W << "\n";  // variable assignment/initialization
    summary << "KE_early=" << KE_early << "\n";  // variable assignment/initialization
    summary << "KE_late=" << KE_late << "\n";  // variable assignment/initialization
    summary << "KE_change_pct=" << KE_change_pct << "\n";  // variable assignment/initialization
    summary << "TE_drift_pct=" << TE_drift_pct << "\n";  // variable assignment/initialization
    summary << "pass_omega=" << pass_omega << "\n";  // variable assignment/initialization
    summary << "pass_heating=" << pass_heating << "\n";  // variable assignment/initialization
    summary << "overall=" << (pass_omega && pass_heating) << "\n";  // variable assignment/initialization
    summary.close();  // statement
// (blank line)
    return (pass_omega && pass_heating) ? 0 : 1;  // return from function
}  // end block
