#include <iostream>  // include iostream
#include <vector>  // include vector
#include <cmath>  // include cmath
#include <fstream>  // include fstream
#include <omp.h>  // include omp.h
#include <algorithm>  // include algorithm
#include <random>  // include random
#include <string>  // include string
#include <mpi.h>  // include mpi.h
// (blank line)
struct Particle { double x, ux, uy, uz, w; int id, species; };  // type definition
constexpr double eps0 = 8.8541878128e-12, me = 9.1093837e-31, qe = -1.602176634e-19, c = 299792458.0;  // variable assignment/initialization
// (blank line)
// Cooley-Tukey Radix-2 Fast Fourier Transform (FFT) Solver
void standalone_fft(std::vector<double>& real, std::vector<double>& imag, int n, bool inverse) {  // begin block
    for (int i = 0, j = 0; i < n; ++i) {  // begin block
        if (i < j) { std::swap(real[i], real[j]); std::swap(imag[i], imag[j]); }  // control flow statement
        int bit = n >> 1;   // variable assignment/initialization
        while (j & bit) { j ^= bit; bit >>= 1; }   // control flow statement
        j ^= bit;  // variable assignment/initialization
    }  // end block
    for (int len = 2; len <= n; len <<= 1) {  // begin block
        double angle = 2.0 * M_PI / len * (inverse ? 1.0 : -1.0);  // function declaration/call or variable initialization
        double wlen_r = std::cos(angle), wlen_i = std::sin(angle);  // function declaration/call or variable initialization
        for (int i = 0; i < n; i += len) {  // begin block
            double w_r = 1.0, w_i = 0.0;  // variable assignment/initialization
            for (int j = 0; j < len / 2; ++j) {  // begin block
                int u = i + j, v = i + j + len / 2;  // variable assignment/initialization
                double t_r = real[v] * w_r - imag[v] * w_i;  // variable assignment/initialization
                double t_i = real[v] * w_i + imag[v] * w_r;  // variable assignment/initialization
                real[v] = real[u] - t_r; imag[v] = imag[u] - t_i;   // variable assignment/initialization
                real[u] += t_r; imag[u] += t_i;  // variable assignment/initialization
                double n_w_r = w_r * wlen_r - w_i * wlen_i;   // variable assignment/initialization
                w_i = w_r * wlen_i + w_i * wlen_r;   // variable assignment/initialization
                w_r = n_w_r;  // variable assignment/initialization
            }  // end block
        }  // end block
    }  // end block
    if (inverse) {   // begin block
        for (int i = 0; i < n; ++i) { real[i] /= n; imag[i] /= n; }   // control flow statement
    }  // end block
}  // end block
// (blank line)
int main(int argc, char* argv[]) {  // begin block
    MPI_Init(&argc, &argv);   // statement
    int mpi_rank, mpi_size;   // statement
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);   // statement
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);  // statement
// (blank line)
    const int Nx_global = 512, ppc_bg = 150, ppc_bm = 50, num_steps = 1500;  // variable assignment/initialization
    const double dx = 2e-5, L_global = Nx_global * dx, dt = 1e-14;  // variable assignment/initialization
    const int Nx_local = Nx_global / mpi_size;  // variable assignment/initialization
    const double x_local_min = mpi_rank * Nx_local * dx, x_local_max = (mpi_rank + 1) * Nx_local * dx;  // function declaration/call or variable initialization
// (blank line)
    std::vector<double> Ex(Nx_local + 4, 0.0), rho(Nx_local + 4, 0.0), local_q(Nx_local + 1, 0.0);  // function declaration/call or variable initialization
    std::vector<Particle> local_particles;  // statement
    const double w_p = (1e19 * L_global) / (Nx_global * ppc_bg), w_b = (1e17 * L_global) / (Nx_global * ppc_bm);  // function declaration/call or variable initialization
// (blank line)
    std::mt19937 gen(42 + mpi_rank);   // function declaration/call or variable initialization
    double Thot = 4000.0, Tcold = 500.0;  // variable assignment/initialization
    double T_start = Thot - (Thot - Tcold) * (static_cast<double>(mpi_rank) / mpi_size);  // function declaration/call or variable initialization
    double T_end   = Thot - (Thot - Tcold) * (static_cast<double>(mpi_rank + 1) / mpi_size);  // function declaration/call or variable initialization
// (blank line)
    for (int i = 0; i < Nx_local * ppc_bg; ++i) {  // begin block
        double Te = T_start * (1.0 - static_cast<double>(i/ppc_bg)/Nx_local) + T_end * (static_cast<double>(i/ppc_bg)/Nx_local);  // function declaration/call or variable initialization
        std::normal_distribution<double> dist_v(0.0, std::sqrt((Te * 1.60218e-19) / me));  // function declaration/call or variable initialization
        Particle p; p.x = x_local_min + ((i / ppc_bg) + 0.5) / ppc_bg * dx; p.w = w_p;  // variable assignment/initialization
        p.ux = dist_v(gen); p.uy = dist_v(gen); p.uz = dist_v(gen); p.species = 0; p.id = i;   // variable assignment/initialization
        local_particles.push_back(p);  // statement
    }  // end block
// (blank line)
    double v_bm = 0.6 * c, g_b = 1.0 / std::sqrt(1.0 - (v_bm * v_bm) / (c * c));  // function declaration/call or variable initialization
    for (int i = 0; i < Nx_local * ppc_bm; ++i) {  // begin block
        Particle p; p.x = x_local_min + ((i / ppc_bm) + 0.5) / ppc_bm * dx; p.w = w_b; p.ux = g_b * v_bm; p.uy = 0.0; p.uz = 0.0; p.species = 1; p.id = i + 100000;  // variable assignment/initialization
        p.x += 0.01 * dx * std::sin(2.0 * M_PI * p.x / L_global);   // variable assignment/initialization
        while (p.x >= L_global) p.x -= L_global;   // control flow statement
        while (p.x < 0.0) p.x += L_global;  // control flow statement
        if (p.x >= x_local_min && p.x < x_local_max) local_particles.push_back(p);  // control flow statement
    }  // end block
// (blank line)
    double* rho_ptr = rho.data();   // function declaration/call or variable initialization
    int left_neighbor = (mpi_rank == 0) ? -1 : mpi_rank - 1;  // function declaration/call or variable initialization
    int right_neighbor = (mpi_rank == mpi_size - 1) ? -1 : mpi_rank + 1;  // function declaration/call or variable initialization
    const int left_proc  = (left_neighbor  == -1) ? MPI_PROC_NULL : left_neighbor;  // function declaration/call or variable initialization
    const int right_proc = (right_neighbor == -1) ? MPI_PROC_NULL : right_neighbor;  // function declaration/call or variable initialization
// (blank line)
    std::ofstream diag_out;   // statement
    if(mpi_rank == 0) {   // begin block
        diag_out.open("gol3_anomalous_conductivity.csv");   // statement
        diag_out << "step,time,max_E,avg_q_parallel\n";   // statement
    }  // end block
// (blank line)
    std::vector<double> rho_global_R(Nx_global, 0.0), rho_global_I(Nx_global, 0.0);  // function declaration/call or variable initialization
    std::vector<double> E_global_R(Nx_global, 0.0), E_global_I(Nx_global, 0.0);  // function declaration/call or variable initialization
// (blank line)
    double t_dep = 0.0, t_comm = 0.0, t_gather = 0.0, t_fft = 0.0,  // code line
           t_scatter = 0.0, t_push = 0.0, t_diag = 0.0, t_total = 0.0;  // variable assignment/initialization
// (blank line)
    for (int step = 1; step <= num_steps; ++step) {  // begin block
        const double t_iter_start = MPI_Wtime();  // function declaration/call or variable initialization
        double t0 = MPI_Wtime();  // function declaration/call or variable initialization
        std::fill(rho.begin(), rho.end(), 0.0);  // function declaration/call or variable initialization
        #pragma omp parallel for reduction(+:rho_ptr[:Nx_local+4]) schedule(static)  // preprocessor directive
        for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
            double x_norm = (local_particles[i].x - x_local_min) / dx;   // function declaration/call or variable initialization
            int j = std::max(0, std::min(static_cast<int>(std::floor(x_norm)), Nx_local - 1));  // function declaration/call or variable initialization
            rho_ptr[j + 1] += (qe * local_particles[i].w * (1.0 - (x_norm - j))) / dx;   // variable assignment/initialization
            rho_ptr[j + 2] += (qe * local_particles[i].w * (x_norm - j)) / dx;  // variable assignment/initialization
        }  // end block
        t_dep += MPI_Wtime() - t0; t0 = MPI_Wtime();  // variable assignment/initialization
// (blank line)
        MPI_Status status;   // statement
        double send_L_ghost = rho[1]; // FIXED: Sends leftmost active node value scalar, not entire vector object pointer 
        double recv_R_ghost = 0.0;  // variable assignment/initialization
// (blank line)
        MPI_Sendrecv(&send_L_ghost, 1, MPI_DOUBLE, left_proc, 0, &recv_R_ghost, 1, MPI_DOUBLE, right_proc, 0, MPI_COMM_WORLD, &status);   // statement
        rho[Nx_local + 1] += recv_R_ghost;   // variable assignment/initialization
// (blank line)
        double local_sum = 0.0, global_sum = 0.0;   // variable assignment/initialization
        for(int j = 1; j <= Nx_local; ++j) local_sum += rho[j];  // control flow statement
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  // statement
        double ion_corr = -global_sum / Nx_global;   // variable assignment/initialization
        for(int j = 1; j <= Nx_local + 1; ++j) rho[j] += ion_corr;  // control flow statement
        t_comm += MPI_Wtime() - t0; t0 = MPI_Wtime();  // variable assignment/initialization
// (blank line)
        std::vector<double> local_rho_pure(Nx_local);   // function declaration/call or variable initialization
        for(int j=0; j<Nx_local; ++j) local_rho_pure[j] = rho[j+1];  // control flow statement
        MPI_Gather(local_rho_pure.data(), Nx_local, MPI_DOUBLE, rho_global_R.data(), Nx_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);  // statement
        t_gather += MPI_Wtime() - t0; t0 = MPI_Wtime();  // variable assignment/initialization
// (blank line)
        if (mpi_rank == 0) {  // begin block
            std::fill(rho_global_I.begin(), rho_global_I.end(), 0.0);   // function declaration/call or variable initialization
            standalone_fft(rho_global_R, rho_global_I, Nx_global, false);  // statement
            std::fill(E_global_R.begin(), E_global_R.end(), 0.0);  // function declaration/call or variable initialization
            std::fill(E_global_I.begin(), E_global_I.end(), 0.0);  // function declaration/call or variable initialization
            #pragma omp parallel for schedule(static)  // preprocessor directive
            for (int k = 1; k < Nx_global; ++k) {   // begin block
                double k_wave = (2.0 / dx) * std::sin(M_PI * k / Nx_global);   // function declaration/call or variable initialization
                E_global_R[k] = -rho_global_I[k] / (eps0 * k_wave);   // variable assignment/initialization
                E_global_I[k] =  rho_global_R[k] / (eps0 * k_wave);   // variable assignment/initialization
            }  // end block
            standalone_fft(E_global_R, E_global_I, Nx_global, true);  // statement
        }  // end block
        t_fft += MPI_Wtime() - t0; t0 = MPI_Wtime();  // variable assignment/initialization
// (blank line)
        std::vector<double> local_E_pure(Nx_local);   // function declaration/call or variable initialization
        MPI_Scatter(E_global_R.data(), Nx_local, MPI_DOUBLE, local_E_pure.data(), Nx_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);  // statement
        for(int j=0; j<Nx_local; ++j) Ex[j+1] = local_E_pure[j];  // control flow statement
        t_scatter += MPI_Wtime() - t0; t0 = MPI_Wtime();  // variable assignment/initialization
// (blank line)
        std::vector<Particle> keep_p, send_L, send_R;  // statement
        #pragma omp parallel  // preprocessor directive
        {  // begin block
            std::vector<Particle> t_keep, t_send_L, t_send_R;  // statement
            #pragma omp for schedule(static)  // preprocessor directive
            for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
                double x_norm = (local_particles[i].x - x_local_min) / dx;   // function declaration/call or variable initialization
                int j = std::max(0, std::min(static_cast<int>(std::floor(x_norm)), Nx_local - 1));  // function declaration/call or variable initialization
                double E = (1.0 - (x_norm - j)) * Ex[j + 1] + (x_norm - j) * Ex[j + 2];  // function declaration/call or variable initialization
                double g = std::sqrt(1.0 + (local_particles[i].ux*local_particles[i].ux + local_particles[i].uy*local_particles[i].uy + local_particles[i].uz*local_particles[i].uz)/(c*c));  // function declaration/call or variable initialization
                local_particles[i].ux -= (qe * dt / me) * E;   // variable assignment/initialization
                local_particles[i].x += (local_particles[i].ux / g) * dt;  // variable assignment/initialization
// (blank line)
                if (local_particles[i].x < x_local_min) t_send_L.push_back(local_particles[i]);  // control flow statement
                else if (local_particles[i].x >= x_local_max) t_send_R.push_back(local_particles[i]);  // else branch
                else t_keep.push_back(local_particles[i]);  // else branch
            }  // end block
            #pragma omp critical  // preprocessor directive
            {   // begin block
                keep_p.insert(keep_p.end(), t_keep.begin(), t_keep.end());   // statement
                send_L.insert(send_L.end(), t_send_L.begin(), t_send_L.end());   // statement
                send_R.insert(send_R.end(), t_send_R.begin(), t_send_R.end());   // statement
            }  // end block
        }  // end block
// (blank line)
        int s_size_L = send_L.size(), r_size_R = 0, s_size_R = send_R.size(), r_size_L = 0;  // function declaration/call or variable initialization
        MPI_Sendrecv(&s_size_L, 1, MPI_INT, left_proc, 1, &r_size_R, 1, MPI_INT, right_proc, 1, MPI_COMM_WORLD, &status);   // statement
        if (left_neighbor == -1) {  // begin block
            for(auto& p : send_L) {   // begin block
                if (p.species == 0) {   // begin block
                    std::normal_distribution<double> d_hot(0.0, std::sqrt((Thot * 1.60218e-19) / me));   // function declaration/call or variable initialization
                    p.x = 0.0 + (static_cast<double>(rand())/RAND_MAX)*dx;   // variable assignment/initialization
                    p.ux = std::abs(d_hot(gen)); p.uy = d_hot(gen); p.uz = d_hot(gen);   // variable assignment/initialization
                    keep_p.push_back(p);   // statement
                }   // end block
            }  // end block
            s_size_L = 0;  // variable assignment/initialization
        }  // end block
// (blank line)
        MPI_Sendrecv(&s_size_R, 1, MPI_INT, right_proc, 2, &r_size_L, 1, MPI_INT, left_proc, 2, MPI_COMM_WORLD, &status);   // statement
        if (right_neighbor == -1) {  // begin block
            for(auto& p : send_R) {   // begin block
                if (p.species == 0) {   // begin block
                    std::normal_distribution<double> d_cold(0.0, std::sqrt((Tcold * 1.60218e-19) / me));   // function declaration/call or variable initialization
                    p.x = L_global - (static_cast<double>(rand())/RAND_MAX)*dx;   // variable assignment/initialization
                    p.ux = -std::abs(d_cold(gen)); p.uy = d_cold(gen); p.uz = d_cold(gen);   // variable assignment/initialization
                    keep_p.push_back(p);   // statement
                }   // end block
            }  // end block
            s_size_R = 0;  // variable assignment/initialization
        }  // end block
// (blank line)
        std::vector<Particle> recv_L(r_size_L), recv_R(r_size_R);  // function declaration/call or variable initialization
        MPI_Sendrecv(send_L.data(), s_size_L*sizeof(Particle), MPI_BYTE, left_proc, 3, recv_R.data(), r_size_R*sizeof(Particle), MPI_BYTE, right_proc, 3, MPI_COMM_WORLD, &status);  // statement
        MPI_Sendrecv(send_R.data(), s_size_R*sizeof(Particle), MPI_BYTE, right_proc, 4, recv_L.data(), r_size_L*sizeof(Particle), MPI_BYTE, left_proc, 4, MPI_COMM_WORLD, &status);  // statement
// (blank line)
        for (auto& p : recv_R) { if (p.x >= L_global) p.x -= L_global; keep_p.push_back(p); }  // control flow statement
        for (auto& p : recv_L) { if (p.x < 0.0) p.x += L_global; keep_p.push_back(p); }  // control flow statement
        local_particles = std::move(keep_p);  // variable assignment/initialization
        t_push += MPI_Wtime() - t0; t0 = MPI_Wtime();  // variable assignment/initialization
// (blank line)
        if (step % 20 == 0) {  // begin block
            std::fill(local_q.begin(), local_q.end(), 0.0);  // function declaration/call or variable initialization
            double* q_ptr = local_q.data();  // function declaration/call or variable initialization
            #pragma omp parallel for reduction(+:q_ptr[:Nx_local+1]) schedule(static)  // preprocessor directive
            for (size_t i = 0; i < local_particles.size(); ++i) {  // begin block
                if (local_particles[i].species == 0) {  // begin block
                    double x_norm = (local_particles[i].x - x_local_min) / dx;  // function declaration/call or variable initialization
                    int j = std::max(0, std::min(static_cast<int>(std::floor(x_norm)), Nx_local - 1));  // function declaration/call or variable initialization
                    double g = std::sqrt(1.0 + (local_particles[i].ux*local_particles[i].ux)/(c*c));  // function declaration/call or variable initialization
                    double vx = local_particles[i].ux / g;  // variable assignment/initialization
                    q_ptr[j] += 0.5 * me * (vx * vx * vx) * local_particles[i].w / dx;  // variable assignment/initialization
                }  // end block
            }  // end block
// (blank line)
            double local_q_sum = 0.0, global_q_avg = 0.0, local_max_E = 0.0, global_max_E = 0.0;  // variable assignment/initialization
            for (int j = 1; j <= Nx_local; ++j) {  // begin block
                local_q_sum += std::abs(local_q[j - 1]);  // variable assignment/initialization
                local_max_E = std::max(local_max_E, std::abs(Ex[j]));  // variable assignment/initialization
            }  // end block
            MPI_Reduce(&local_q_sum, &global_q_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  // statement
            MPI_Reduce(&local_max_E, &global_max_E, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);  // statement
// (blank line)
            if (mpi_rank == 0) {  // begin block
                diag_out << step << "," << step*dt << "," << global_max_E << "," << (global_q_avg / Nx_global) << "\n";  // statement
            }  // end block
        }  // end block
        t_diag += MPI_Wtime() - t0;  // variable assignment/initialization
        t_total += MPI_Wtime() - t_iter_start;  // variable assignment/initialization
// (blank line)
        if (step % 20 == 0 || step == 1) {  // begin block
            long long local_n = static_cast<long long>(local_particles.size()), global_n = 0;  // function declaration/call or variable initialization
            MPI_Reduce(&local_n, &global_n, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);  // statement
            double t_total_max = 0.0;  // variable assignment/initialization
            MPI_Reduce(&t_total, &t_total_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);  // statement
            if (mpi_rank == 0) {  // begin block
                double per = t_total / step;  // variable assignment/initialization
                std::cout << "[step " << step << "/" << num_steps << "] "  // code line
                          << "avg=" << (1e3 * per) << " ms/step  ETA=" << (per * (num_steps - step)) << " s  "  // code line
                          << "particles(global)=" << global_n << "  max_wall_rank=" << t_total_max << " s\n"  // code line
                          << "    breakdown ms/step: dep=" << (1e3 * t_dep / step)  // code line
                          << " comm=" << (1e3 * t_comm / step)  // code line
                          << " gather=" << (1e3 * t_gather / step)  // code line
                          << " fft=" << (1e3 * t_fft / step)  // code line
                          << " scatter=" << (1e3 * t_scatter / step)  // code line
                          << " push=" << (1e3 * t_push / step)  // code line
                          << " diag=" << (1e3 * t_diag / step) << std::endl;  // variable assignment/initialization
            }  // end block
        }  // end block
    }  // end block
    if (mpi_rank == 0) {  // begin block
        diag_out.close();  // statement
        std::cout << "\n=== phase totals over " << num_steps << " steps (rank 0, seconds) ===\n"  // code line
                  << "  total (loop) : " << t_total << "\n"  // code line
                  << "  deposit      : " << t_dep << "\n"  // code line
                  << "  comm (ghost/neutrality): " << t_comm << "\n"  // code line
                  << "  gather       : " << t_gather << "\n"  // code line
                  << "  fft (rank 0) : " << t_fft << "\n"  // code line
                  << "  scatter      : " << t_scatter << "\n"  // code line
                  << "  push+move    : " << t_push << "\n"  // code line
                  << "  diagnostics  : " << t_diag << "\n";  // statement
    }  // end block
    if (mpi_rank == 0) {  // begin block
        diag_out.close();  // statement
        std::cout << "\n====================================================\n";  // variable assignment/initialization
        std::cout << " GOL-3 NON-LINEAR SCALING RUN COMPLETION            \n";  // statement
        std::cout << "====================================================\n";  // variable assignment/initialization
        std::cout << " Status: Complete. Telemetry saved to gol3_anomalous_conductivity.csv\n";  // statement
        std::cout << "====================================================\n";  // variable assignment/initialization
    }  // end block
    MPI_Finalize();  // statement
    return 0;  // return from function
}  // end block
